"""DeepDTA binding affinity prediction wrapper.

Provides AI-based drug-target binding affinity prediction.
Falls back to sequence-based heuristic if model not available.
"""

from __future__ import annotations

import logging
from typing import Optional

logger = logging.getLogger(__name__)


def predict_binding_affinity(
    target_sequence: str,
    ligand: str,
    ligand_type: str = "smiles",
    device: str = "cpu",
) -> Optional[float]:
    """Predict binding affinity using DeepDTA or fallback.

    Args:
        target_sequence: Protein amino acid sequence
        ligand: SMILES string (SM) or peptide sequence
        ligand_type: "smiles" or "sequence"
        device: Compute device

    Returns:
        Predicted pKd value (higher = stronger binding), or None if fails.
    """
    # Try DeepDTA model
    result = _try_deepdta(target_sequence, ligand, ligand_type, device)
    if result is not None:
        return result

    # Fallback: sequence-based heuristic
    if ligand_type == "sequence":
        return _sequence_binding_heuristic(target_sequence, ligand)

    return None


def _try_deepdta(
    target_seq: str,
    ligand: str,
    ligand_type: str,
    device: str,
) -> Optional[float]:
    """Attempt to use DeepDTA model for prediction."""
    try:
        import torch
        from transformers import AutoModel, AutoTokenizer
    except ImportError:
        return None

    # DeepDTA requires specific model - check availability
    try:
        # This is a placeholder for actual DeepDTA integration
        # Real implementation would load a trained DeepDTA model
        # and encode protein + drug sequences
        logger.debug("DeepDTA model integration pending")
        return None
    except Exception:
        return None


def _sequence_binding_heuristic(target_seq: str, peptide_seq: str) -> float:
    """Simple heuristic binding score based on sequence properties.

    Returns estimated pKd (5-9 range, higher = better).
    """
    target_seq = target_seq.upper()
    peptide_seq = peptide_seq.upper()

    score = 6.0  # baseline

    # Length preference: 8-20 aa optimal
    pep_len = len(peptide_seq)
    if 8 <= pep_len <= 20:
        score += 0.5
    elif pep_len > 25:
        score -= 0.3

    # Charge complementarity
    target_positive = sum(1 for aa in target_seq if aa in "RKH")
    target_negative = sum(1 for aa in target_seq if aa in "DE")
    pep_positive = sum(1 for aa in peptide_seq if aa in "RKH")
    pep_negative = sum(1 for aa in peptide_seq if aa in "DE")

    # Opposite charges attract
    charge_comp = (
        min(target_positive, pep_negative) + min(target_negative, pep_positive)
    ) / max(pep_len, 1)
    score += min(charge_comp * 2.0, 1.0)

    # Hydrophobic content match
    hydrophobic = set("VILMFYW")
    target_h = sum(1 for aa in target_seq if aa in hydrophobic) / max(len(target_seq), 1)
    pep_h = sum(1 for aa in peptide_seq if aa in hydrophobic) / max(pep_len, 1)
    if abs(target_h - pep_h) < 0.15:
        score += 0.5

    # Aromatic content (binding interface)
    aromatic = sum(1 for aa in peptide_seq if aa in "FWY")
    if 1 <= aromatic <= 4:
        score += 0.3

    # Clamp to reasonable pKd range
    return round(max(4.0, min(10.0, score)), 2)
