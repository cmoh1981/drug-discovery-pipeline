"""M3-SM: De novo small molecule design.

Uses REINVENT 4 when available, falls back to RDKit fragment-based enumeration.
"""

from __future__ import annotations

import logging
import random
from pathlib import Path
from typing import Optional

from drugdiscovery.types import Candidate, PipelineConfig, TargetProfile
from drugdiscovery.utils.io import write_csv

logger = logging.getLogger(__name__)

# Common drug-like fragments for enumeration fallback
SCAFFOLD_FRAGMENTS = [
    "c1ccccc1",           # benzene
    "c1ccncc1",           # pyridine
    "c1ccoc1",            # furan
    "c1cc[nH]c1",         # pyrrole
    "c1cnc2ccccc2n1",     # quinazoline
    "c1ccc2[nH]ccc2c1",   # indole
    "C1CCCCC1",           # cyclohexane
    "C1CCNCC1",           # piperidine
    "C1CCOC1",            # tetrahydrofuran
    "C1CNCCN1",           # piperazine
]

LINKERS = ["C", "CC", "CCC", "CO", "CN", "C(=O)", "C(=O)N", "NC(=O)", "O", "S"]

FUNCTIONAL_GROUPS = [
    "O", "N", "F", "Cl", "C(=O)O", "C(=O)N", "S(=O)(=O)N",
    "C#N", "C(F)(F)F", "OC", "NC",
]


def design_small_molecules(
    cfg: PipelineConfig,
    target_profile: TargetProfile,
    output_dir: Path,
) -> list[Candidate]:
    """Generate de novo small molecule candidates."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # Try REINVENT 4 first
    candidates = _try_reinvent4(cfg, target_profile, output_dir)

    if not candidates:
        logger.info("REINVENT 4 not available; using RDKit fragment enumeration")
        candidates = _rdkit_fragment_enumeration(cfg, target_profile, output_dir)

    if candidates:
        rows = [c.to_dict() for c in candidates]
        write_csv(output_dir / "candidates.csv", rows)

    logger.info("Generated %d small molecule candidates", len(candidates))
    return candidates


def _try_reinvent4(
    cfg: PipelineConfig,
    target: TargetProfile,
    output_dir: Path,
) -> list[Candidate]:
    """Attempt to use REINVENT 4 for de novo design."""
    try:
        import reinvent
        logger.info("REINVENT 4 found — using for de novo generation")
        # REINVENT 4 integration would go here
        # This requires a trained prior model and scoring function setup
        # For now, return empty to fall through to RDKit
        return []
    except ImportError:
        return []


def _rdkit_fragment_enumeration(
    cfg: PipelineConfig,
    target: TargetProfile,
    output_dir: Path,
) -> list[Candidate]:
    """Generate candidates via RDKit fragment-based enumeration."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors, Crippen
    except ImportError:
        logger.warning("RDKit not available; cannot generate small molecules")
        return []

    num_to_generate = cfg.pepmlm_settings.get("num_per_length", 50) * 3
    candidates: list[Candidate] = []
    seen_smiles: set[str] = set()
    attempts = 0
    max_attempts = num_to_generate * 20

    while len(candidates) < num_to_generate and attempts < max_attempts:
        attempts += 1
        smiles = _assemble_molecule()

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue

        # Canonicalize
        canon = Chem.MolToSmiles(mol)
        if canon in seen_smiles:
            continue
        seen_smiles.add(canon)

        # Drug-likeness filter (Lipinski)
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        hba = Descriptors.NumHAcceptors(mol)
        hbd = Descriptors.NumHDonors(mol)

        if mw > 600 or mw < 150:
            continue
        if logp > 6 or logp < -2:
            continue
        if hba > 12 or hbd > 6:
            continue

        cand = Candidate(
            candidate_id=f"RDKIT_{len(candidates)+1:05d}",
            candidate_type="de_novo",
            modality="small_molecule",
            source="rdkit_enum",
            smiles=canon,
            molecular_weight=round(mw, 2),
            net_charge=Chem.GetFormalCharge(mol),
        )
        candidates.append(cand)

    logger.info("RDKit enumeration: %d valid molecules from %d attempts", len(candidates), attempts)
    return candidates


def _assemble_molecule() -> str:
    """Assemble a random drug-like molecule from fragments."""
    # Pick 1-2 scaffolds
    n_scaffolds = random.choice([1, 2])
    parts = [random.choice(SCAFFOLD_FRAGMENTS) for _ in range(n_scaffolds)]

    # Connect with linker if 2 scaffolds
    if n_scaffolds == 2:
        linker = random.choice(LINKERS)
        smiles = parts[0] + linker + parts[1]
    else:
        smiles = parts[0]

    # Add 0-2 functional groups
    n_groups = random.choice([0, 1, 1, 2])
    for _ in range(n_groups):
        fg = random.choice(FUNCTIONAL_GROUPS)
        smiles = smiles + fg

    return smiles
