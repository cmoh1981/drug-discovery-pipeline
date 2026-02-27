"""Amino acid property tables and physicochemical calculators.

Refactored from YARS2 pipeline: score_candidates.py, phase5_optimize.py, phase6_admet.py
"""

from __future__ import annotations

import math
import re
from typing import Optional

# ---------------------------------------------------------------------------
# Amino acid property tables
# ---------------------------------------------------------------------------

AA_1TO3 = {
    "A": "Ala", "C": "Cys", "D": "Asp", "E": "Glu", "F": "Phe",
    "G": "Gly", "H": "His", "I": "Ile", "K": "Lys", "L": "Leu",
    "M": "Met", "N": "Asn", "P": "Pro", "Q": "Gln", "R": "Arg",
    "S": "Ser", "T": "Thr", "V": "Val", "W": "Trp", "Y": "Tyr",
}

STANDARD_AA = set(AA_1TO3.keys())

AA_MW: dict[str, float] = {
    "A": 89.09, "R": 174.20, "N": 132.12, "D": 133.10, "C": 121.16,
    "E": 147.13, "Q": 146.15, "G": 75.03, "H": 155.16, "I": 131.17,
    "L": 131.17, "K": 146.19, "M": 149.21, "F": 165.19, "P": 115.13,
    "S": 105.09, "T": 119.12, "V": 117.15, "W": 204.23, "Y": 181.19,
}

AA_HYDROPHOBICITY: dict[str, float] = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "E": -3.5, "Q": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
}

# pKa values: (pKa_amino, pKa_carboxyl, pKa_sidechain_or_None)
AA_PKA: dict[str, tuple[float, float, Optional[float]]] = {
    "A": (9.87, 2.35, None), "R": (9.09, 2.18, 12.48),
    "N": (8.80, 2.02, None), "D": (9.82, 1.88, 3.65),
    "C": (10.78, 1.71, 8.33), "E": (9.67, 2.19, 4.25),
    "Q": (9.13, 2.17, None), "G": (9.60, 2.34, None),
    "H": (9.17, 1.82, 6.00), "I": (9.76, 2.32, None),
    "L": (9.60, 2.36, None), "K": (8.95, 2.20, 10.53),
    "M": (9.21, 2.28, None), "F": (9.24, 2.58, None),
    "P": (10.60, 1.99, None), "S": (9.15, 2.21, None),
    "T": (9.12, 2.09, None), "V": (9.72, 2.29, None),
    "W": (9.39, 2.38, None), "Y": (9.11, 2.20, 10.07),
}

# Charged residues at pH 7.4
POSITIVE_RESIDUES = {"R", "K", "H"}
NEGATIVE_RESIDUES = {"D", "E"}
HELIX_PRONE = {"A", "E", "L", "K", "M", "Q", "R"}

# Protease cleavage site patterns
PROTEASE_SITES = {
    "trypsin": re.compile(r"[KR](?!P)"),
    "chymotrypsin": re.compile(r"[FWY](?!P)"),
    "glu_c": re.compile(r"[DE]"),
}

# Instability index dipeptide weights (Guruprasad et al. 1990, representative subset)
DIWV: dict[str, float] = {
    "WW": 1.0, "WC": 1.0, "WM": 24.68, "WH": 24.68,
    "DG": 1.0, "DP": -6.54, "DW": 1.0, "DA": 1.0,
    "EG": 1.0, "EA": 1.0, "ED": 1.0, "EE": 1.0,
    "GG": 1.0, "GD": 1.0, "GN": -7.49, "GA": 1.0,
    "KK": 1.0, "KD": 1.0, "KR": 1.0, "KG": 1.0,
    "RK": 1.0, "RR": 1.0, "RD": 1.0, "RG": 1.0,
    "LL": 1.0, "LI": 1.0, "LV": 1.0, "LA": 1.0,
    "II": 1.0, "IV": 1.0, "IL": 1.0, "IA": 1.0,
    "VV": 1.0, "VL": 1.0, "VI": 1.0, "VA": 1.0,
    "FF": 1.0, "FW": 1.0, "FY": 1.0, "FA": 1.0,
    "AA": 1.0, "AL": 1.0, "AV": 1.0, "AI": 1.0,
}


# ---------------------------------------------------------------------------
# Physicochemical calculators
# ---------------------------------------------------------------------------

def compute_molecular_weight(sequence: str) -> float:
    """Compute molecular weight in Daltons from amino acid sequence."""
    seq = _clean_sequence(sequence)
    if not seq:
        return 0.0
    mw = sum(AA_MW.get(aa, 110.0) for aa in seq)
    # Subtract water for each peptide bond
    mw -= 18.015 * (len(seq) - 1)
    return round(mw, 2)


def compute_net_charge(sequence: str, ph: float = 7.4) -> float:
    """Compute net charge at given pH using Henderson-Hasselbalch."""
    seq = _clean_sequence(sequence)
    if not seq:
        return 0.0

    charge = 0.0
    # N-terminal amino group (use first residue pKa)
    pka_nh2 = AA_PKA.get(seq[0], (9.0, 2.0, None))[0]
    charge += 1.0 / (1.0 + 10 ** (ph - pka_nh2))

    # C-terminal carboxyl group
    pka_cooh = AA_PKA.get(seq[-1], (9.0, 2.0, None))[1]
    charge -= 1.0 / (1.0 + 10 ** (pka_cooh - ph))

    # Side chains
    for aa in seq:
        pka_data = AA_PKA.get(aa)
        if pka_data and pka_data[2] is not None:
            pka_sc = pka_data[2]
            if aa in POSITIVE_RESIDUES:
                charge += 1.0 / (1.0 + 10 ** (ph - pka_sc))
            elif aa in NEGATIVE_RESIDUES or aa in ("C", "Y"):
                charge -= 1.0 / (1.0 + 10 ** (pka_sc - ph))

    return round(charge, 2)


def compute_gravy(sequence: str) -> float:
    """Grand Average of Hydropathicity (Kyte-Doolittle)."""
    seq = _clean_sequence(sequence)
    if not seq:
        return 0.0
    total = sum(AA_HYDROPHOBICITY.get(aa, 0.0) for aa in seq)
    return round(total / len(seq), 3)


def compute_isoelectric_point(sequence: str, precision: float = 0.01) -> float:
    """Compute pI by bisection on net_charge."""
    lo, hi = 0.0, 14.0
    for _ in range(100):
        mid = (lo + hi) / 2.0
        charge = compute_net_charge(sequence, ph=mid)
        if abs(charge) < precision:
            return round(mid, 2)
        if charge > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2.0, 2)


def compute_instability_index(sequence: str) -> float:
    """Instability index (Guruprasad et al. 1990). >40 = unstable."""
    seq = _clean_sequence(sequence)
    if len(seq) < 2:
        return 0.0
    total = 0.0
    for i in range(len(seq) - 1):
        dipeptide = seq[i] + seq[i + 1]
        total += DIWV.get(dipeptide, 1.0)
    return round((10.0 / len(seq)) * total, 2)


def compute_aromaticity(sequence: str) -> float:
    """Fraction of aromatic residues (F, W, Y)."""
    seq = _clean_sequence(sequence)
    if not seq:
        return 0.0
    aromatic = sum(1 for aa in seq if aa in "FWY")
    return round(aromatic / len(seq), 4)


def count_protease_sites(sequence: str) -> dict[str, int]:
    """Count predicted protease cleavage sites."""
    seq = _clean_sequence(sequence)
    return {name: len(pat.findall(seq)) for name, pat in PROTEASE_SITES.items()}


def check_drug_likeness(sequence: str) -> tuple[float, list[str]]:
    """Five-point drug-likeness filter for peptides.

    Returns (score 0-100, list of violations).
    """
    seq = _clean_sequence(sequence)
    violations = []

    # 1. Length: 8-25 aa
    if not (8 <= len(seq) <= 25):
        violations.append(f"length={len(seq)} (want 8-25)")

    # 2. Net charge: [-2, +4]
    charge = compute_net_charge(seq)
    if not (-2 <= charge <= 4):
        violations.append(f"charge={charge:.1f} (want -2 to +4)")

    # 3. No aggregation motifs: 4+ consecutive hydrophobic
    if re.search(r"[VILMFYW]{4,}", seq):
        violations.append("hydrophobic stretch >= 4 residues")

    # 4. Protease sites <= 3
    sites = count_protease_sites(seq)
    total_sites = sites.get("trypsin", 0) + sites.get("chymotrypsin", 0)
    if total_sites > 3:
        violations.append(f"protease_sites={total_sites} (want <=3)")

    # 5. No homopolymers > 3
    if re.search(r"(.)\1{3,}", seq):
        violations.append("homopolymer > 3 consecutive identical residues")

    score = max(0.0, 100.0 - 20.0 * len(violations))
    return score, violations


def compute_peptide_properties(sequence: str) -> dict:
    """Compute all physicochemical properties for a peptide sequence."""
    seq = _clean_sequence(sequence)
    dl_score, dl_violations = check_drug_likeness(seq)
    return {
        "sequence": seq,
        "length": len(seq),
        "molecular_weight": compute_molecular_weight(seq),
        "net_charge": compute_net_charge(seq),
        "gravy": compute_gravy(seq),
        "isoelectric_point": compute_isoelectric_point(seq),
        "instability_index": compute_instability_index(seq),
        "aromaticity": compute_aromaticity(seq),
        "drug_likeness": dl_score,
        "drug_likeness_violations": dl_violations,
        "protease_sites": count_protease_sites(seq),
    }


def _clean_sequence(sequence: str) -> str:
    """Normalize a peptide sequence: uppercase, strip, standard AA only."""
    seq = sequence.upper().strip()
    return "".join(c for c in seq if c in STANDARD_AA)


def safe_float(value, default: float = 0.0) -> float:
    """Safely convert a value to float."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return default
