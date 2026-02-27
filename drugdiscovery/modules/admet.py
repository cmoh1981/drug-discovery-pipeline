"""M7: ADMET prediction module.

Supports both peptide-specific local computation and small molecule
ADMET via ADMETlab 3.0 REST API.

Refactored from YARS2 pipeline: phase6_admet.py
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Optional

from drugdiscovery.types import ADMETProfile, Candidate, PipelineConfig
from drugdiscovery.utils.chemistry import (
    STANDARD_AA,
    compute_aromaticity,
    compute_gravy,
    compute_instability_index,
    compute_molecular_weight,
    compute_net_charge,
    safe_float,
)
from drugdiscovery.utils.io import write_csv

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def predict_admet(
    cfg: PipelineConfig,
    candidates: list[Candidate],
    output_dir: Path,
) -> list[Candidate]:
    """Predict ADMET profiles for all candidates."""
    output_dir.mkdir(parents=True, exist_ok=True)
    profiles: list[dict] = []

    for cand in candidates:
        if cand.modality == "peptide":
            profile = _predict_peptide_admet(cand)
        else:
            profile = _predict_sm_admet(cand)

        cand.admet_score = profile.aggregate_score
        cand.admet_flags = profile.flag_count
        profiles.append(_profile_to_dict(profile))

    write_csv(output_dir / "admet_profiles.csv", profiles)
    logger.info("ADMET profiles saved for %d candidates", len(profiles))
    return candidates


# ---------------------------------------------------------------------------
# Peptide ADMET (local computation)
# ---------------------------------------------------------------------------

def _predict_peptide_admet(cand: Candidate) -> ADMETProfile:
    """Compute ADMET properties for a peptide candidate."""
    seq = cand.sequence
    profile = ADMETProfile(candidate_id=cand.candidate_id)

    # Solubility (GRAVY-based)
    gravy = compute_gravy(seq)
    profile.solubility = _gravy_to_solubility(gravy)

    # Cell-penetrating peptide score
    profile.cpp_score = _predict_cpp_score(seq)
    profile.permeability = profile.cpp_score

    # Protease stability
    profile.protease_stability = _predict_protease_stability(seq)

    # Hemolysis risk
    profile.hemolysis_risk = _predict_hemolysis_risk(seq)

    # Aggregation propensity
    profile.aggregation_propensity = _predict_aggregation(seq)

    # BBB permeability
    profile.bbb_permeability = _predict_bbb(seq)

    # Half-life estimate
    mw = compute_molecular_weight(seq)
    profile.half_life_estimate = "short" if mw < 1500 else "moderate"

    # Compute flags
    flags = []
    if profile.solubility < 0.3:
        flags.append("low_solubility")
    if profile.cpp_score < 0.3:
        flags.append("poor_cell_penetration")
    if profile.protease_stability < 0.4:
        flags.append("protease_susceptible")
    if profile.hemolysis_risk > 0.6:
        flags.append("hemolysis_risk")
    if profile.aggregation_propensity > 0.6:
        flags.append("aggregation_prone")
    if profile.bbb_permeability > 0.5:
        flags.append("bbb_permeable")  # May be wanted or unwanted

    profile.flags = flags
    profile.flag_count = len([f for f in flags if f != "bbb_permeable"])
    profile.aggregate_score = _compute_aggregate(profile)

    return profile


def _gravy_to_solubility(gravy: float) -> float:
    """Convert GRAVY to solubility score (0-1). Negative GRAVY = hydrophilic = good."""
    if gravy <= -2.0:
        return 1.0
    if gravy >= 2.0:
        return 0.0
    return round(0.5 - gravy / 4.0, 4)


def _predict_cpp_score(seq: str) -> float:
    """Predict cell-penetrating peptide likelihood."""
    seq = seq.upper()
    length = len(seq)
    if length == 0:
        return 0.0

    score = 0.0
    # Arginine content (strong CPP indicator)
    r_frac = seq.count("R") / length
    score += min(r_frac * 3.0, 0.4)

    # Positive charge
    charge = compute_net_charge(seq)
    if 2 <= charge <= 8:
        score += 0.2
    elif charge > 0:
        score += 0.1

    # Amphipathicity (alternating hydrophobic/hydrophilic)
    hydrophobic = set("VILMFYW")
    alternations = sum(
        1 for i in range(length - 1)
        if (seq[i] in hydrophobic) != (seq[i + 1] in hydrophobic)
    )
    amphipathic_ratio = alternations / max(length - 1, 1)
    score += amphipathic_ratio * 0.2

    # Length penalty (CPPs typically 5-30 aa)
    if 5 <= length <= 30:
        score += 0.1
    elif length > 30:
        score -= 0.1

    # Tryptophan content (membrane interaction)
    w_frac = seq.count("W") / length
    score += min(w_frac * 2.0, 0.1)

    return round(max(0.0, min(1.0, score)), 4)


def _predict_protease_stability(seq: str) -> float:
    """Predict protease stability (1 = very stable, 0 = unstable)."""
    seq = seq.upper()
    length = len(seq)
    if length == 0:
        return 0.0

    # Count cleavage sites
    trypsin_sites = len(re.findall(r"[KR](?!P)", seq))
    chymo_sites = len(re.findall(r"[FWY](?!P)", seq))
    total_sites = trypsin_sites + chymo_sites

    # Normalize: fewer sites = more stable
    site_density = total_sites / length
    stability = max(0.0, 1.0 - site_density * 3.0)

    # Check for modifications (lowercase = D-amino acid or N-methyl)
    has_mods = any(c.islower() for c in seq)
    if has_mods:
        stability = min(1.0, stability + 0.2)

    return round(stability, 4)


def _predict_hemolysis_risk(seq: str) -> float:
    """Predict hemolysis risk (0 = safe, 1 = high risk)."""
    seq = seq.upper()
    length = len(seq)
    if length == 0:
        return 0.0

    risk = 0.0

    # High positive charge + hydrophobicity = hemolytic
    charge = compute_net_charge(seq)
    gravy = compute_gravy(seq)

    if charge > 4 and gravy > 0:
        risk += 0.4
    elif charge > 2 and gravy > 0.5:
        risk += 0.3

    # Amphipathic helical potential
    hydrophobic = set("VILMFYW")
    h_frac = sum(1 for aa in seq if aa in hydrophobic) / length
    if h_frac > 0.4 and charge > 2:
        risk += 0.3

    # Length: longer peptides more hemolytic
    if length > 20:
        risk += 0.1

    return round(min(1.0, risk), 4)


def _predict_aggregation(seq: str) -> float:
    """Predict aggregation propensity."""
    seq = seq.upper()
    if not seq:
        return 0.0

    score = 0.0

    # Hydrophobic stretches
    if re.search(r"[VILMFYW]{4,}", seq):
        score += 0.4
    if re.search(r"[VILMFYW]{6,}", seq):
        score += 0.3

    # Overall hydrophobicity
    gravy = compute_gravy(seq)
    if gravy > 1.0:
        score += 0.2
    elif gravy > 0.5:
        score += 0.1

    # Homopolymers
    if re.search(r"(.)\1{3,}", seq):
        score += 0.2

    return round(min(1.0, score), 4)


def _predict_bbb(seq: str) -> float:
    """Predict blood-brain barrier permeability."""
    seq = seq.upper()
    if not seq:
        return 0.0

    mw = compute_molecular_weight(seq)
    charge = compute_net_charge(seq)
    gravy = compute_gravy(seq)

    score = 0.0
    # Small, lipophilic, low charge favors BBB crossing
    if mw < 500:
        score += 0.3
    elif mw < 1000:
        score += 0.1
    if -1 <= charge <= 1:
        score += 0.2
    if gravy > 0:
        score += 0.2
    if len(seq) <= 10:
        score += 0.1

    return round(min(1.0, score), 4)


# ---------------------------------------------------------------------------
# Small Molecule ADMET (ADMETlab 3.0 API)
# ---------------------------------------------------------------------------

def _predict_sm_admet(cand: Candidate) -> ADMETProfile:
    """Compute ADMET for small molecules.

    Tries ADMETlab 3.0 API first, falls back to RDKit descriptors.
    """
    profile = ADMETProfile(candidate_id=cand.candidate_id)

    # Try ADMETlab 3.0 API first
    if cand.smiles:
        try:
            from drugdiscovery.tools.admetlab3 import predict_admet_admetlab3

            api_result = predict_admet_admetlab3(cand.smiles)
            if api_result is not None:
                profile = _admetlab3_to_profile(cand.candidate_id, api_result)
                logger.info("[M7] ADMETlab 3.0 prediction for %s", cand.candidate_id)
                return profile
        except Exception as exc:
            logger.warning("[M7] ADMETlab 3.0 failed for %s: %s", cand.candidate_id, exc)

    # Fall back to RDKit-based local computation
    try:
        profile = _rdkit_admet(cand.smiles, profile)
    except Exception as exc:
        logger.warning("RDKit ADMET failed for %s: %s", cand.candidate_id, exc)
        profile.aggregate_score = 0.5  # neutral fallback

    return profile


def _rdkit_admet(smiles: str, profile: ADMETProfile) -> ADMETProfile:
    """Compute ADMET properties using RDKit descriptors."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Crippen, Lipinski
    except ImportError:
        logger.warning("RDKit not available for SM ADMET")
        profile.aggregate_score = 0.5
        return profile

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.warning("Invalid SMILES: %s", smiles)
        profile.aggregate_score = 0.3
        return profile

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hba = Descriptors.NumHAcceptors(mol)
    hbd = Descriptors.NumHDonors(mol)
    tpsa = Descriptors.TPSA(mol)
    rotatable = Descriptors.NumRotatableBonds(mol)

    flags = []

    # Lipinski Rule of 5
    lipinski_violations = 0
    if mw > 500:
        lipinski_violations += 1
    if logp > 5:
        lipinski_violations += 1
    if hba > 10:
        lipinski_violations += 1
    if hbd > 5:
        lipinski_violations += 1

    # Oral bioavailability
    profile.oral_bioavailability = max(0.0, 1.0 - lipinski_violations * 0.25)

    # Solubility (logP-based)
    if logp < 0:
        profile.solubility = 0.9
    elif logp < 3:
        profile.solubility = 0.7
    elif logp < 5:
        profile.solubility = 0.4
    else:
        profile.solubility = 0.1
        flags.append("low_solubility")

    # Permeability (TPSA-based)
    if tpsa < 90:
        profile.permeability = 0.8
    elif tpsa < 140:
        profile.permeability = 0.5
    else:
        profile.permeability = 0.2
        flags.append("poor_permeability")

    # BBB (TPSA < 60 and logP 1-3)
    if tpsa < 60 and 1 < logp < 3:
        profile.bbb_permeability = 0.7
    else:
        profile.bbb_permeability = 0.2

    # hERG liability (crude: highly lipophilic + basic amine)
    if logp > 3.5:
        profile.herg_liability = 0.4
        flags.append("herg_risk")
    else:
        profile.herg_liability = 0.1

    # Hepatotoxicity (MW > 600 or logP > 4)
    if mw > 600 or logp > 4:
        profile.hepatotoxicity_risk = 0.4
        flags.append("hepatotoxicity_risk")
    else:
        profile.hepatotoxicity_risk = 0.1

    profile.flags = flags
    profile.flag_count = len(flags)
    profile.aggregate_score = _compute_aggregate(profile)

    return profile


# ---------------------------------------------------------------------------
# ADMETlab 3.0 result conversion
# ---------------------------------------------------------------------------

def _clamp_float(val, lo: float = 0.0, hi: float = 1.0) -> float:
    """Clamp a value to [lo, hi], handling None and non-numeric values."""
    try:
        v = float(val) if val is not None else 0.0
    except (TypeError, ValueError):
        v = 0.0
    return max(lo, min(hi, v))


def _admetlab3_to_profile(candidate_id: str, data: dict) -> ADMETProfile:
    """Convert ADMETlab 3.0 API results to ADMETProfile."""
    profile = ADMETProfile(candidate_id=candidate_id)

    # Solubility: logS → 0-1 score (higher logS = more soluble = better)
    logs = data.get("solubility", 0.0)
    if logs >= 0:
        profile.solubility = 1.0
    elif logs >= -2:
        profile.solubility = 0.8
    elif logs >= -4:
        profile.solubility = 0.5
    elif logs >= -6:
        profile.solubility = 0.2
    else:
        profile.solubility = 0.0

    # Permeability from Caco-2 (log Papp; > -5 is good)
    caco2 = data.get("caco2_permeability", 0.0)
    profile.permeability = _clamp_float(0.5 + caco2 / 10.0)

    # Oral bioavailability (classification probability)
    profile.oral_bioavailability = _clamp_float(data.get("oral_bioavailability", 0.5))

    # BBB penetration (classification probability)
    profile.bbb_permeability = _clamp_float(data.get("bbb_penetration", 0.0))

    # hERG blocker probability (high = risky)
    profile.herg_liability = _clamp_float(data.get("herg_blocker", 0.0))

    # Hepatotoxicity / DILI probability (high = risky)
    profile.hepatotoxicity_risk = _clamp_float(data.get("hepatotoxicity", 0.0))

    # Build flags
    flags: list[str] = []
    if profile.solubility < 0.3:
        flags.append("low_solubility")
    if profile.permeability < 0.3:
        flags.append("poor_permeability")
    if profile.herg_liability > 0.5:
        flags.append("herg_risk")
    if profile.hepatotoxicity_risk > 0.5:
        flags.append("hepatotoxicity_risk")
    if _clamp_float(data.get("ames_toxicity", 0.0)) > 0.5:
        flags.append("ames_positive")
    if _clamp_float(data.get("carcinogenicity", 0.0)) > 0.5:
        flags.append("carcinogenicity_risk")

    profile.flags = flags
    profile.flag_count = len(flags)
    profile.aggregate_score = _compute_aggregate(profile)

    return profile


# ---------------------------------------------------------------------------
# Aggregate scoring
# ---------------------------------------------------------------------------

def _compute_aggregate(profile: ADMETProfile) -> float:
    """Compute weighted aggregate ADMET score (0-1, higher = better)."""
    components = [
        (profile.solubility, 0.20),
        (profile.permeability, 0.15),
        (1.0 - profile.hemolysis_risk, 0.15),
        (profile.protease_stability, 0.15),
        (1.0 - profile.aggregation_propensity, 0.10),
        (1.0 - profile.hepatotoxicity_risk, 0.10),
        (1.0 - profile.herg_liability, 0.10),
        (profile.oral_bioavailability, 0.05),
    ]
    score = sum(val * wt for val, wt in components)
    return round(max(0.0, min(1.0, score)), 4)


def _profile_to_dict(profile: ADMETProfile) -> dict:
    """Convert ADMETProfile to flat dict for CSV."""
    return {
        "candidate_id": profile.candidate_id,
        "solubility": round(profile.solubility, 4),
        "permeability": round(profile.permeability, 4),
        "cpp_score": round(profile.cpp_score, 4),
        "oral_bioavailability": round(profile.oral_bioavailability, 4),
        "bbb_permeability": round(profile.bbb_permeability, 4),
        "protease_stability": round(profile.protease_stability, 4),
        "hemolysis_risk": round(profile.hemolysis_risk, 4),
        "aggregation_propensity": round(profile.aggregation_propensity, 4),
        "hepatotoxicity_risk": round(profile.hepatotoxicity_risk, 4),
        "herg_liability": round(profile.herg_liability, 4),
        "flag_count": profile.flag_count,
        "flags": "; ".join(profile.flags),
        "aggregate_score": round(profile.aggregate_score, 4),
    }
