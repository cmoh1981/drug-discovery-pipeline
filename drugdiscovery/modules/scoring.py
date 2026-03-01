"""M5 – Binding & Scoring module.

Computes physicochemical properties, drug-likeness, selectivity, sequence
diversity, and a weighted composite score for every candidate.  Saves
scored_candidates.csv and top_N_candidates.csv to *output_dir*.
"""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import Any

from drugdiscovery.types import PipelineConfig, Candidate, TargetProfile
from drugdiscovery.utils.chemistry import (
    compute_molecular_weight,
    compute_net_charge,
    compute_gravy,
    compute_isoelectric_point,
    compute_instability_index,
    check_drug_likeness,
    safe_float,
    STANDARD_AA,
)
from drugdiscovery.utils.io import write_csv

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _levenshtein_normalized(s1: str, s2: str) -> float:
    """Normalized Levenshtein distance in [0, 1] (1 = maximally different).

    Returns 1.0 if either sequence is empty (treat as fully different).
    """
    if not s1 or not s2:
        return 1.0
    if s1 == s2:
        return 0.0

    len1, len2 = len(s1), len(s2)

    # Allocate two rows for DP to keep memory O(min(m,n))
    if len1 < len2:
        s1, s2 = s2, s1
        len1, len2 = len2, len1

    # s1 is the longer string
    prev = list(range(len2 + 1))
    curr = [0] * (len2 + 1)

    for i in range(1, len1 + 1):
        curr[0] = i
        for j in range(1, len2 + 1):
            cost = 0 if s1[i - 1] == s2[j - 1] else 1
            curr[j] = min(
                curr[j - 1] + 1,        # insertion
                prev[j] + 1,            # deletion
                prev[j - 1] + cost,     # substitution
            )
        prev, curr = curr, prev

    distance = prev[len2]
    return round(distance / max(len1, len2), 4)


def _composition_similarity(seq1: str, seq2: str) -> float:
    """Cosine similarity of amino-acid composition vectors.

    Each vector has one dimension per standard amino acid (20-dim), holding
    the fractional count of that amino acid in the sequence.  Returns a
    value in [0, 1]; returns 0.0 if either sequence is empty.
    """
    if not seq1 or not seq2:
        return 0.0

    aa_list = sorted(STANDARD_AA)  # fixed ordering

    def _vec(seq: str) -> list[float]:
        n = len(seq) or 1
        return [seq.count(aa) / n for aa in aa_list]

    v1 = _vec(seq1)
    v2 = _vec(seq2)

    dot = sum(a * b for a, b in zip(v1, v2))
    norm1 = math.sqrt(sum(a * a for a in v1))
    norm2 = math.sqrt(sum(b * b for b in v2))

    if norm1 == 0.0 or norm2 == 0.0:
        return 0.0

    return round(max(0.0, min(1.0, dot / (norm1 * norm2))), 4)


def _pocket_sequence(pocket_residue_names: list[str]) -> str:
    """Convert a list of three-letter residue names to a one-letter sequence.

    Unknown residue names are skipped.
    """
    three_to_one = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    }
    result = []
    for name in pocket_residue_names:
        upper = name.upper().strip()
        if upper in three_to_one:
            result.append(three_to_one[upper])
        elif len(upper) == 1 and upper in STANDARD_AA:
            # Already single-letter
            result.append(upper)
    return "".join(result)


def _clamp(value: float, lo: float = 0.0, hi: float = 1.0) -> float:
    return max(lo, min(hi, value))


# ---------------------------------------------------------------------------
# Step 1 – Physicochemical scoring
# ---------------------------------------------------------------------------

def score_physicochemical(candidate: Candidate) -> Candidate:
    """Populate MW, net_charge, gravy, isoelectric_point on *candidate*.

    Peptides: computed analytically from sequence.
    Small molecules: attempted via RDKit from SMILES; graceful fallback.
    """
    if candidate.modality == "peptide":
        seq = candidate.sequence
        if not seq:
            logger.warning(
                "[M5] Candidate %s has no sequence; skipping physicochemical scoring.",
                candidate.candidate_id,
            )
            return candidate

        candidate.molecular_weight = safe_float(compute_molecular_weight(seq))
        candidate.net_charge = safe_float(compute_net_charge(seq))
        candidate.gravy = safe_float(compute_gravy(seq))
        candidate.isoelectric_point = safe_float(compute_isoelectric_point(seq))

        logger.debug(
            "[M5] %s physicochemical: MW=%.2f, charge=%.2f, GRAVY=%.3f, pI=%.2f",
            candidate.candidate_id,
            candidate.molecular_weight,
            candidate.net_charge,
            candidate.gravy,
            candidate.isoelectric_point,
        )

    else:  # small molecule
        smiles = candidate.smiles
        if not smiles:
            logger.warning(
                "[M5] Candidate %s has no SMILES; skipping physicochemical scoring.",
                candidate.candidate_id,
            )
            return candidate

        try:
            from rdkit import Chem  # type: ignore
            from rdkit.Chem import Descriptors  # type: ignore

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logger.warning(
                    "[M5] RDKit could not parse SMILES for %s; skipping SM physicochemical.",
                    candidate.candidate_id,
                )
                return candidate

            candidate.molecular_weight = round(Descriptors.ExactMolWt(mol), 4)
            # net_charge approximated from formal charges
            candidate.net_charge = round(
                safe_float(Chem.GetFormalCharge(mol)), 4
            )
            # GRAVY / pI are peptide-specific; leave at defaults for SM
            logger.debug(
                "[M5] %s SM physicochemical: MW=%.2f, charge=%.2f",
                candidate.candidate_id,
                candidate.molecular_weight,
                candidate.net_charge,
            )

        except ImportError:
            logger.warning(
                "[M5] RDKit not available; SM physicochemical scoring skipped for %s.",
                candidate.candidate_id,
            )

    return candidate


# ---------------------------------------------------------------------------
# Step 2 – Drug-likeness scoring
# ---------------------------------------------------------------------------

def score_drug_likeness(candidate: Candidate) -> Candidate:
    """Assign drug_likeness score (0-1) to *candidate*.

    Peptides: use check_drug_likeness(); normalize 0-100 → 0-1.
    Small molecules: Lipinski Rule of 5 via RDKit; fallback 0.5.
    """
    if candidate.modality == "peptide":
        seq = candidate.sequence
        if not seq:
            candidate.drug_likeness = 0.0
            return candidate

        raw_score, violations = check_drug_likeness(seq)
        candidate.drug_likeness = round(_clamp(raw_score / 100.0), 4)

        if violations:
            logger.debug(
                "[M5] %s drug-likeness violations: %s",
                candidate.candidate_id,
                "; ".join(violations),
            )

    else:  # small molecule – Lipinski Ro5
        smiles = candidate.smiles
        if not smiles:
            candidate.drug_likeness = 0.5
            return candidate

        try:
            from rdkit import Chem  # type: ignore
            from rdkit.Chem import Descriptors, rdMolDescriptors  # type: ignore

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logger.warning(
                    "[M5] RDKit parse failed for %s; drug_likeness defaulted to 0.5.",
                    candidate.candidate_id,
                )
                candidate.drug_likeness = 0.5
                return candidate

            mw = Descriptors.ExactMolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = rdMolDescriptors.CalcNumHBD(mol)
            hba = rdMolDescriptors.CalcNumHBA(mol)

            violations = sum([
                mw > 500,
                logp > 5,
                hbd > 5,
                hba > 10,
            ])
            # Each violation costs 0.25 of the score; 0 violations → 1.0
            candidate.drug_likeness = round(_clamp(1.0 - 0.25 * violations), 4)

            logger.debug(
                "[M5] %s Ro5: MW=%.1f logP=%.2f HBD=%d HBA=%d → drug_likeness=%.4f",
                candidate.candidate_id,
                mw, logp, hbd, hba,
                candidate.drug_likeness,
            )

        except ImportError:
            logger.warning(
                "[M5] RDKit not available; drug_likeness defaulted to 0.5 for %s.",
                candidate.candidate_id,
            )
            candidate.drug_likeness = 0.5

    return candidate


# ---------------------------------------------------------------------------
# Step 3 – Selectivity scoring
# ---------------------------------------------------------------------------

def score_selectivity(
    candidate: Candidate,
    target_profile: TargetProfile,
) -> Candidate:
    """Assign selectivity_score (0-1) to *candidate*.

    Peptides:
      * Build a pseudo-sequence from target binding pocket residue names.
      * Compute composition-vector cosine similarity vs. that pocket.
      * If anti-target binding pockets exist, compute similarity vs. those.
      * selectivity = sim_target / (sim_target + sim_anti_target)

    Small molecules: placeholder 0.5 until structure-based docking is wired in.
    """
    if candidate.modality == "peptide":
        seq = candidate.sequence
        if not seq:
            candidate.selectivity_score = 0.0
            return candidate

        # Build target pocket sequence from first binding pocket that has residue names
        target_pocket_seq = ""
        for pocket in target_profile.binding_pockets:
            if pocket.residue_names:
                target_pocket_seq = _pocket_sequence(pocket.residue_names)
                break
        # Fallback: use a substring of the full protein sequence centred on known sites
        if not target_pocket_seq and target_profile.sequence:
            target_pocket_seq = target_profile.sequence[:50]

        sim_target = _composition_similarity(seq, target_pocket_seq) if target_pocket_seq else 0.5

        # Build anti-target pocket sequence(s)
        anti_sims: list[float] = []
        for anti in target_profile.anti_targets:
            # anti_targets are free-form dicts; try common key names for pocket residues
            pocket_residues: list[str] = []
            for key in ("residue_names", "binding_residues", "residues", "pocket_residues"):
                if key in anti and isinstance(anti[key], list):
                    pocket_residues = anti[key]
                    break

            anti_seq = _pocket_sequence(pocket_residues) if pocket_residues else ""
            if not anti_seq:
                # Fallback to sequence field if available
                anti_seq = anti.get("sequence", "")[:50] if isinstance(anti.get("sequence"), str) else ""

            if anti_seq:
                anti_sims.append(_composition_similarity(seq, anti_seq))

        if anti_sims:
            sim_anti = sum(anti_sims) / len(anti_sims)
        else:
            # No anti-target info → assume moderate anti-target similarity
            sim_anti = 0.3

        denom = sim_target + sim_anti
        if denom < 1e-9:
            selectivity = 0.5
        else:
            selectivity = sim_target / denom

        candidate.selectivity_score = round(_clamp(selectivity), 4)

        logger.debug(
            "[M5] %s selectivity: sim_target=%.4f sim_anti=%.4f → %.4f",
            candidate.candidate_id,
            sim_target,
            sim_anti,
            candidate.selectivity_score,
        )

    else:  # small molecule – placeholder
        candidate.selectivity_score = 0.5
        logger.debug(
            "[M5] %s (SM) selectivity set to placeholder 0.5", candidate.candidate_id
        )

    return candidate


# ---------------------------------------------------------------------------
# Step 4 – Sequence diversity
# ---------------------------------------------------------------------------

def compute_sequence_diversity(candidates: list[Candidate]) -> list[Candidate]:
    """Assign sequence_diversity (0-1) to every candidate in the batch.

    For peptides: mean pairwise normalized Levenshtein distance against all
    other peptide candidates.  A candidate that is very different from all
    others scores close to 1.0; one identical to most others scores close to
    0.0.

    For small molecules: Tanimoto diversity using RDKit Morgan fingerprints
    (radius 2).  Falls back to 0.5 if RDKit is unavailable.
    """
    peptides = [c for c in candidates if c.modality == "peptide"]
    sms = [c for c in candidates if c.modality != "peptide"]

    # ---- Peptide diversity ----
    MAX_DIVERSITY_CANDIDATES = 200
    if len(peptides) > 1:
        sequences = [c.sequence for c in peptides]
        n = len(sequences)
        if n > MAX_DIVERSITY_CANDIDATES:
            import random
            logger.info(
                "[M5] Sampling %d/%d candidates for diversity (O(n^2) cap)",
                MAX_DIVERSITY_CANDIDATES, n,
            )
        for i, cand in enumerate(peptides):
            if n > MAX_DIVERSITY_CANDIDATES:
                sample_indices = random.sample(
                    [j for j in range(n) if j != i],
                    min(MAX_DIVERSITY_CANDIDATES, n - 1),
                )
            else:
                sample_indices = [j for j in range(n) if j != i]
            distances = [
                _levenshtein_normalized(sequences[i], sequences[j])
                for j in sample_indices
            ]
            mean_dist = sum(distances) / len(distances)
            cand.sequence_diversity = round(_clamp(mean_dist), 4)
            logger.debug(
                "[M5] %s sequence_diversity=%.4f (peptide, n=%d)",
                cand.candidate_id,
                cand.sequence_diversity,
                n,
            )
    elif len(peptides) == 1:
        # Single peptide – maximum diversity by default (no competitors)
        peptides[0].sequence_diversity = 1.0

    # ---- Small-molecule diversity ----
    if sms:
        try:
            from rdkit import Chem  # type: ignore
            from rdkit.Chem import rdMolDescriptors, DataStructs  # type: ignore

            fps: list[Any] = []
            valid_sms: list[Candidate] = []
            for cand in sms:
                mol = Chem.MolFromSmiles(cand.smiles) if cand.smiles else None
                if mol is not None:
                    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
                    fps.append(fp)
                    valid_sms.append(cand)
                else:
                    cand.sequence_diversity = 0.5

            n = len(valid_sms)
            if n > 1:
                for i, cand in enumerate(valid_sms):
                    similarities = [
                        DataStructs.TanimotoSimilarity(fps[i], fps[j])
                        for j in range(n)
                        if j != i
                    ]
                    mean_sim = sum(similarities) / len(similarities)
                    # Diversity = 1 - mean similarity
                    cand.sequence_diversity = round(_clamp(1.0 - mean_sim), 4)
                    logger.debug(
                        "[M5] %s sequence_diversity=%.4f (SM Tanimoto)",
                        cand.candidate_id,
                        cand.sequence_diversity,
                    )
            elif n == 1:
                valid_sms[0].sequence_diversity = 1.0

        except ImportError:
            logger.warning(
                "[M5] RDKit not available; SM sequence_diversity defaulted to 0.5 for %d candidates.",
                len(sms),
            )
            for cand in sms:
                cand.sequence_diversity = 0.5

    return candidates


# ---------------------------------------------------------------------------
# Step 5 – Composite scoring
# ---------------------------------------------------------------------------

# Absolute binding energy threshold (kcal/mol).  Candidates weaker than this
# receive zero binding credit.  More negative = stronger binding.
BINDING_THRESHOLD = -5.0

# Hard discard threshold (kcal/mol).  Candidates with binding_score above
# this value are removed entirely before ranking — they are too weak to
# be considered drug candidates.  Only applied when real docking data exists
# (i.e. at least some candidates have binding_score < 0).
BINDING_DISCARD_THRESHOLD = -3.0


def compute_composite_score(
    candidate: Candidate,
    weights: dict,
    *,
    binding_min: float,
    binding_max: float,
    target_mode: str = "antagonist",
) -> Candidate:
    """Compute weighted composite score (0-1) for a single candidate.

    Scored components (weighted):
      binding_energy   – lower raw binding_score is better; min-max normalize
                         then invert so that the best binder → 1.0.
                         Candidates weaker than BINDING_THRESHOLD get 0.0.
      selectivity      – already 0-1 (higher is better)
      drug_likeness    – already 0-1 (higher is better)
      admet_aggregate  – already 0-1 (higher is better)
      moa_consistency  – 1.0 if predicted MOA matches target, 0.5 unknown

    Gates (not scored, applied as multipliers):
      structure_confidence – if pLDDT/iPTM < 0.5, composite is halved

    Post-ranking filters (not in composite):
      sequence_diversity – applied separately after ranking

    *binding_min* / *binding_max* are the batch-level min/max of raw
    binding_score values; they must be pre-computed and passed in.
    """
    # ---- Absolute binding threshold ----
    if candidate.binding_score > BINDING_THRESHOLD:
        norm_binding = 0.0
    else:
        # ---- Normalize binding score (batch-relative) ----
        binding_range = binding_max - binding_min
        if binding_range < 1e-9:
            norm_binding = 0.5
        else:
            raw_norm = (candidate.binding_score - binding_min) / binding_range
            norm_binding = 1.0 - raw_norm
        norm_binding = _clamp(norm_binding)

    # ---- MOA consistency ----
    moa_raw = candidate.moa_predicted
    if not isinstance(moa_raw, str):
        moa_raw = "unknown"
    moa = (moa_raw or "unknown").lower()
    if moa == target_mode.lower():
        moa_score = 1.0
    elif moa == "unknown":
        moa_score = 0.5
    else:
        moa_score = 0.0

    components = {
        "binding_energy": norm_binding,
        "selectivity": _clamp(safe_float(candidate.selectivity_score)),
        "drug_likeness": _clamp(safe_float(candidate.drug_likeness)),
        "admet_aggregate": _clamp(safe_float(candidate.admet_score)),
        "moa_consistency": moa_score,
        "perturbation": _clamp(safe_float(candidate.perturbation_score)),
    }

    total_weight = sum(weights.get(k, 0.0) for k in components)
    if total_weight < 1e-9:
        logger.warning(
            "[M5] Scoring weights sum to zero for %s; composite_score=0.0",
            candidate.candidate_id,
        )
        candidate.composite_score = 0.0
        return candidate

    composite = sum(
        weights.get(k, 0.0) * v for k, v in components.items()
    ) / total_weight

    # ---- Structure confidence gate ----
    conf = safe_float(candidate.structure_confidence)
    if 0 < conf < 0.5:
        composite *= 0.5
        logger.debug(
            "[M5] %s structure_confidence=%.2f < 0.5; composite halved",
            candidate.candidate_id, conf,
        )

    candidate.composite_score = round(_clamp(composite), 4)

    logger.debug(
        "[M5] %s composite=%.4f "
        "(bind=%.4f sel=%.4f dl=%.4f admet=%.4f moa=%.4f pert=%.4f conf_gate=%.2f)",
        candidate.candidate_id,
        candidate.composite_score,
        components["binding_energy"],
        components["selectivity"],
        components["drug_likeness"],
        components["admet_aggregate"],
        components["moa_consistency"],
        components["perturbation"],
        conf,
    )

    return candidate


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def compute_sa_score(candidate: Candidate) -> Candidate:
    """Compute Synthetic Accessibility score using RDKit (1-10, lower = easier).

    Only applies to small molecules with valid SMILES.
    """
    if candidate.modality == "peptide" or not candidate.smiles:
        return candidate
    try:
        from rdkit import Chem
        from rdkit.Chem import RDConfig
        import sys
        import os

        sa_path = os.path.join(RDConfig.RDContribDir, "SA_Score")
        if sa_path not in sys.path:
            sys.path.insert(0, sa_path)
        import sascorer  # type: ignore

        mol = Chem.MolFromSmiles(candidate.smiles)
        if mol is not None:
            candidate.sa_score = round(sascorer.calculateScore(mol), 2)
            logger.debug(
                "[M5] %s SA score: %.2f",
                candidate.candidate_id,
                candidate.sa_score,
            )
    except Exception as exc:
        logger.debug("[M5] SA score computation failed for %s: %s", candidate.candidate_id, exc)

    return candidate


def _apply_binding_discard_filter(candidates: list[Candidate]) -> list[Candidate]:
    """Remove candidates with binding_score above BINDING_DISCARD_THRESHOLD.

    Only applies when real docking data is present (at least some candidates
    have non-zero binding scores).  If all binding scores are 0.0 (no docking
    was performed), all candidates are kept.
    """
    has_real_docking = any(c.binding_score != 0.0 for c in candidates)
    if not has_real_docking:
        logger.info("[M5] No real docking data; skipping binding discard filter")
        return candidates

    before = len(candidates)
    kept = [c for c in candidates if c.binding_score <= BINDING_DISCARD_THRESHOLD]
    discarded = before - len(kept)

    if discarded > 0:
        logger.info(
            "[M5] Binding discard filter: removed %d/%d candidates "
            "(binding_score > %.1f kcal/mol)",
            discarded, before, BINDING_DISCARD_THRESHOLD,
        )

    # Safety: never discard ALL candidates
    if not kept:
        logger.warning(
            "[M5] All candidates would be discarded by binding filter; "
            "keeping top 10 by binding_score"
        )
        kept = sorted(candidates, key=lambda c: c.binding_score)[:10]

    return kept


def _apply_diversity_filter(candidates: list[Candidate], top_n: int) -> list[Candidate]:
    """Post-ranking diversity filter: remove near-duplicates from top-N.

    For small molecules: if two candidates have Tanimoto similarity > 0.85,
    keep the one with the higher composite score.
    For peptides: if two candidates have normalized Levenshtein distance < 0.15,
    keep the higher-scoring one.

    Applied AFTER ranking, only to the top-N selection.
    """
    if len(candidates) <= 1:
        return candidates

    top = candidates[:top_n]
    rest = candidates[top_n:]

    sm_top = [c for c in top if c.modality != "peptide"]
    pep_top = [c for c in top if c.modality == "peptide"]

    # --- SM diversity filter via Tanimoto ---
    if len(sm_top) > 1:
        try:
            from rdkit import Chem
            from rdkit.Chem import rdMolDescriptors, DataStructs

            fps = []
            valid = []
            for c in sm_top:
                mol = Chem.MolFromSmiles(c.smiles) if c.smiles else None
                if mol is not None:
                    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
                    fps.append(fp)
                    valid.append(c)
                else:
                    valid.append(c)
                    fps.append(None)

            to_remove: set[int] = set()
            for i in range(len(valid)):
                if i in to_remove or fps[i] is None:
                    continue
                for j in range(i + 1, len(valid)):
                    if j in to_remove or fps[j] is None:
                        continue
                    sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
                    if sim > 0.85:
                        # Remove the lower-scoring one (later in list = lower score)
                        to_remove.add(j)

            if to_remove:
                sm_top = [c for idx, c in enumerate(valid) if idx not in to_remove]
                logger.info(
                    "[M5] Diversity filter removed %d near-duplicate SM candidates",
                    len(to_remove),
                )
        except ImportError:
            pass

    # --- Peptide diversity filter via Levenshtein ---
    if len(pep_top) > 1:
        to_remove_pep: set[int] = set()
        for i in range(len(pep_top)):
            if i in to_remove_pep:
                continue
            for j in range(i + 1, len(pep_top)):
                if j in to_remove_pep:
                    continue
                dist = _levenshtein_normalized(pep_top[i].sequence, pep_top[j].sequence)
                if dist < 0.15:
                    to_remove_pep.add(j)
        if to_remove_pep:
            pep_top = [c for idx, c in enumerate(pep_top) if idx not in to_remove_pep]
            logger.info(
                "[M5] Diversity filter removed %d near-duplicate peptide candidates",
                len(to_remove_pep),
            )

    # Recombine and backfill from rest if we removed candidates
    diverse_top = sorted(sm_top + pep_top, key=lambda c: c.composite_score, reverse=True)

    # Backfill from rest to maintain top_n count
    while len(diverse_top) < top_n and rest:
        diverse_top.append(rest.pop(0))

    return diverse_top + rest


def score_candidates(
    cfg: PipelineConfig,
    candidates: list[Candidate],
    target_profile: TargetProfile,
    output_dir: Path,
) -> list[Candidate]:
    """Score all candidates and save results to *output_dir*.

    Steps performed (per candidate):
      1. score_physicochemical  – MW, charge, GRAVY, pI
      2. score_drug_likeness    – drug_likeness ∈ [0,1]
      3. score_selectivity      – selectivity_score ∈ [0,1]
      4. compute_sa_score       – SA score (1-10) for SM

    Then for the full batch:
      5. _apply_binding_discard_filter – remove extremely weak binders
      6. compute_sequence_diversity    – sequence_diversity ∈ [0,1]
      7. compute_composite_score       – composite_score ∈ [0,1]

    Finally:
      8. Sort by composite_score descending, assign rank
      9. _apply_diversity_filter – remove near-duplicates from top-N
     10. Write scored_candidates.csv and top_N_candidates.csv

    Returns the sorted, scored candidate list.
    """
    if not candidates:
        logger.warning("[M5] score_candidates called with empty candidate list.")
        return candidates

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("[M5] Scoring %d candidates …", len(candidates))

    # --- Per-candidate passes ---
    for i, cand in enumerate(candidates):
        logger.debug("[M5] Processing candidate %d/%d: %s", i + 1, len(candidates), cand.candidate_id)
        score_physicochemical(cand)
        score_drug_likeness(cand)
        score_selectivity(cand, target_profile)
        compute_sa_score(cand)

    # --- Binding discard filter (remove extremely weak binders) ---
    candidates = _apply_binding_discard_filter(candidates)

    # --- Batch: sequence diversity ---
    logger.info("[M5] Computing sequence diversity …")
    candidates = compute_sequence_diversity(candidates)

    # --- Batch: composite score ---
    weights = cfg.scoring_weights

    binding_scores = [c.binding_score for c in candidates]
    binding_min = min(binding_scores)
    binding_max = max(binding_scores)

    logger.info(
        "[M5] Binding score range: [%.4f, %.4f]",
        binding_min,
        binding_max,
    )

    for cand in candidates:
        compute_composite_score(
            cand,
            weights,
            binding_min=binding_min,
            binding_max=binding_max,
            target_mode=cfg.mode.value,
        )

    # --- Sort and rank ---
    candidates.sort(key=lambda c: c.composite_score, reverse=True)

    # --- Post-ranking diversity filter ---
    top_n = cfg.top_n
    candidates = _apply_diversity_filter(candidates, top_n)

    for rank, cand in enumerate(candidates, start=1):
        cand.rank = rank

    # --- Save outputs ---
    all_rows = [c.to_dict() for c in candidates]
    write_csv(output_dir / "scored_candidates.csv", all_rows)
    logger.info("[M5] Saved scored_candidates.csv (%d rows)", len(all_rows))

    top_rows = all_rows[:top_n]
    write_csv(output_dir / f"top_{top_n}_candidates.csv", top_rows)
    logger.info("[M5] Saved top_%d_candidates.csv (%d rows)", top_n, len(top_rows))

    logger.info(
        "[M5] Scoring complete. Best: %s (composite=%.4f)",
        candidates[0].candidate_id,
        candidates[0].composite_score,
    )

    return candidates
