"""M6: Candidate modification strategies.

Peptide: cyclization, stapling, D-amino acid, N-methylation, MPP, PEGylation
Small molecule: prodrug, bioisostere, salt form

Refactored from YARS2 pipeline: phase5_optimize.py
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

from drugdiscovery.types import Candidate, PipelineConfig
from drugdiscovery.utils.chemistry import (
    HELIX_PRONE,
    STANDARD_AA,
    compute_gravy,
    compute_isoelectric_point,
    compute_instability_index,
    compute_molecular_weight,
    compute_net_charge,
    check_drug_likeness,
)
from drugdiscovery.utils.io import write_csv

logger = logging.getLogger(__name__)

# MPP motifs for mitochondrial targeting
MPP_MOTIFS = {
    "SS31": "RDFK",
    "FxR3": "FRFRFR",
    "MTS_short": "RRFF",
}


def modify_candidates(
    cfg: PipelineConfig,
    candidates: list[Candidate],
    output_dir: Path,
) -> list[Candidate]:
    """Apply modification strategies and return new variant candidates."""
    output_dir.mkdir(parents=True, exist_ok=True)
    modified: list[Candidate] = []

    if cfg.modality.value == "peptide":
        strategies = cfg.pepmlm_settings.get("modifications", [
            "cyclization", "stapling", "d_amino_acid",
            "n_methylation", "mpp_conjugation", "pegylation",
        ])
        for cand in candidates:
            variants = _apply_peptide_modifications(cand, strategies)
            modified.extend(variants)
    else:
        for cand in candidates:
            variants = _apply_sm_modifications(cand)
            modified.extend(variants)

    # Save modification log
    if modified:
        rows = [_variant_to_dict(v) for v in modified]
        write_csv(output_dir / "modified_candidates.csv", rows)

    logger.info("Generated %d modified variants from %d parents", len(modified), len(candidates))
    return modified


# ---------------------------------------------------------------------------
# Peptide modifications
# ---------------------------------------------------------------------------

def _apply_peptide_modifications(
    parent: Candidate,
    strategies: list[str],
) -> list[Candidate]:
    """Apply all applicable peptide modifications to a parent candidate."""
    variants: list[Candidate] = []
    seq = parent.sequence.upper()
    counter = 0

    if "cyclization" in strategies:
        for v in _cyclization(parent, seq):
            counter += 1
            v.candidate_id = f"{parent.candidate_id}_CYC{counter}"
            variants.append(v)

    if "stapling" in strategies:
        for v in _stapling(parent, seq):
            counter += 1
            v.candidate_id = f"{parent.candidate_id}_STP{counter}"
            variants.append(v)

    if "d_amino_acid" in strategies:
        for v in _d_amino_acid(parent, seq):
            counter += 1
            v.candidate_id = f"{parent.candidate_id}_DAA{counter}"
            variants.append(v)

    if "n_methylation" in strategies:
        for v in _n_methylation(parent, seq):
            counter += 1
            v.candidate_id = f"{parent.candidate_id}_NME{counter}"
            variants.append(v)

    if "mpp_conjugation" in strategies:
        for v in _mpp_conjugation(parent, seq):
            counter += 1
            v.candidate_id = f"{parent.candidate_id}_MPP{counter}"
            variants.append(v)

    if "pegylation" in strategies:
        for v in _pegylation(parent, seq):
            counter += 1
            v.candidate_id = f"{parent.candidate_id}_PEG{counter}"
            variants.append(v)

    return variants


def _make_variant(parent: Candidate, seq: str, mod_type: str, detail: str) -> Candidate:
    """Create a modified variant from a parent candidate."""
    dl_score, _ = check_drug_likeness(seq)
    return Candidate(
        candidate_id="",  # Set by caller
        candidate_type="modified",
        modality=parent.modality,
        source=parent.source,
        sequence=seq,
        smiles=parent.smiles,
        molecular_weight=compute_molecular_weight(seq),
        net_charge=compute_net_charge(seq),
        gravy=compute_gravy(seq),
        isoelectric_point=compute_isoelectric_point(seq),
        drug_likeness=dl_score / 100.0,
        parent_id=parent.candidate_id,
        modification=mod_type,
        modification_detail=detail,
    )


def _cyclization(parent: Candidate, seq: str) -> list[Candidate]:
    """Head-to-tail and disulfide cyclization."""
    variants = []

    # Head-to-tail: feasible if length >= 8 and no Pro at positions 0-1
    if len(seq) >= 8 and seq[0] != "P" and seq[1] != "P":
        detail = "Head-to-tail cyclization (amide bond N-to-C terminus)"
        variants.append(_make_variant(parent, seq, "cyclization", detail))

    # Disulfide: add Cys at both termini
    if "C" not in seq[:2] and "C" not in seq[-2:]:
        ds_seq = "C" + seq + "C"
        detail = "Disulfide cyclization (Cys added at both termini)"
        variants.append(_make_variant(parent, ds_seq, "cyclization_disulfide", detail))

    return variants


def _stapling(parent: Candidate, seq: str) -> list[Candidate]:
    """Hydrocarbon stapling at i, i+4 helical positions."""
    variants = []

    # Find helical runs >= 7 residues
    run_start = None
    run_len = 0
    for i, aa in enumerate(seq):
        if aa in HELIX_PRONE:
            if run_start is None:
                run_start = i
            run_len = i - run_start + 1
        else:
            if run_len >= 7 and run_start is not None:
                # Staple at run_start and run_start + 4
                pos_i = run_start
                pos_j = run_start + 4
                detail = f"Hydrocarbon staple at positions {pos_i+1} and {pos_j+1} (i, i+4)"
                variants.append(_make_variant(parent, seq, "stapling", detail))
            run_start = None
            run_len = 0

    # Check final run
    if run_len >= 7 and run_start is not None:
        pos_i = run_start
        pos_j = run_start + 4
        detail = f"Hydrocarbon staple at positions {pos_i+1} and {pos_j+1} (i, i+4)"
        variants.append(_make_variant(parent, seq, "stapling", detail))

    return variants


def _d_amino_acid(parent: Candidate, seq: str) -> list[Candidate]:
    """D-amino acid substitution at protease cleavage sites."""
    variants = []

    # Find trypsin sites (K/R not before P)
    for match in re.finditer(r"[KR](?!P)", seq):
        pos = match.start()
        new_seq = seq[:pos] + seq[pos].lower() + seq[pos + 1:]
        detail = f"D-{seq[pos]} at position {pos+1} (trypsin resistance)"
        variants.append(_make_variant(parent, new_seq, "d_amino_acid", detail))

    # Find chymotrypsin sites (F/W/Y not before P) - limit to first 2
    chymo_count = 0
    for match in re.finditer(r"[FWY](?!P)", seq):
        if chymo_count >= 2:
            break
        pos = match.start()
        new_seq = seq[:pos] + seq[pos].lower() + seq[pos + 1:]
        detail = f"D-{seq[pos]} at position {pos+1} (chymotrypsin resistance)"
        variants.append(_make_variant(parent, new_seq, "d_amino_acid", detail))
        chymo_count += 1

    return variants


def _n_methylation(parent: Candidate, seq: str) -> list[Candidate]:
    """N-methylation at non-helical positions."""
    # Find non-helical positions
    candidates_pos = [i for i, aa in enumerate(seq) if aa not in HELIX_PRONE and aa != "P"]

    if not candidates_pos:
        return []

    # Select up to 3 evenly spaced positions
    if len(candidates_pos) <= 3:
        selected = candidates_pos
    else:
        step = len(candidates_pos) / 3
        selected = [candidates_pos[int(i * step)] for i in range(3)]

    new_seq = list(seq)
    positions = []
    for pos in selected:
        new_seq[pos] = new_seq[pos].lower()
        positions.append(str(pos + 1))

    mod_seq = "".join(new_seq)
    detail = f"N-methylation at positions {', '.join(positions)}"
    return [_make_variant(parent, mod_seq, "n_methylation", detail)]


def _mpp_conjugation(parent: Candidate, seq: str) -> list[Candidate]:
    """Mitochondria-penetrating peptide conjugation."""
    variants = []
    for name, motif in MPP_MOTIFS.items():
        conj_seq = seq + motif
        detail = f"{name} motif ({motif}) appended to C-terminus"
        variants.append(_make_variant(parent, conj_seq, "mpp_conjugation", detail))
    return variants


def _pegylation(parent: Candidate, seq: str) -> list[Candidate]:
    """PEGylation strategy annotation.

    Since PEG is not an amino acid, we annotate the candidate
    and adjust MW estimate rather than modifying the sequence.
    """
    # Find Lys or Cys for PEG attachment
    peg_sites = []
    for i, aa in enumerate(seq):
        if aa in ("K", "C"):
            peg_sites.append((i, aa))

    if not peg_sites:
        return []

    # Use first suitable site
    site_pos, site_aa = peg_sites[0]
    peg_mw_2k = 2000.0
    peg_mw_5k = 5000.0

    variants = []
    for peg_size, peg_mw in [("PEG2k", peg_mw_2k), ("PEG5k", peg_mw_5k)]:
        v = _make_variant(parent, seq, "pegylation", f"{peg_size} at {site_aa}{site_pos+1}")
        v.molecular_weight += peg_mw
        v.metadata["peg_size"] = peg_size
        v.metadata["peg_site"] = f"{site_aa}{site_pos+1}"
        variants.append(v)

    return variants


# ---------------------------------------------------------------------------
# Small molecule modifications
# ---------------------------------------------------------------------------

def _apply_sm_modifications(parent: Candidate) -> list[Candidate]:
    """Apply small molecule modification strategies (annotation-based)."""
    variants = []

    # Prodrug annotation
    v = Candidate(
        candidate_id=f"{parent.candidate_id}_PRODRUG1",
        candidate_type="modified",
        modality=parent.modality,
        source=parent.source,
        smiles=parent.smiles,
        molecular_weight=parent.molecular_weight,
        net_charge=parent.net_charge,
        parent_id=parent.candidate_id,
        modification="prodrug",
        modification_detail="Ester prodrug strategy (requires synthesis)",
    )
    variants.append(v)

    # Salt form annotation
    v2 = Candidate(
        candidate_id=f"{parent.candidate_id}_SALT1",
        candidate_type="modified",
        modality=parent.modality,
        source=parent.source,
        smiles=parent.smiles,
        molecular_weight=parent.molecular_weight,
        net_charge=parent.net_charge,
        parent_id=parent.candidate_id,
        modification="salt_form",
        modification_detail="HCl salt form for improved solubility",
    )
    variants.append(v2)

    return variants


def _variant_to_dict(v: Candidate) -> dict:
    """Convert variant to dict for CSV export."""
    d = v.to_dict()
    d["parent_id"] = v.parent_id
    d["modification"] = v.modification
    d["modification_detail"] = v.modification_detail
    return d
