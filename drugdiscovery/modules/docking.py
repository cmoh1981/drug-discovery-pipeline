"""M4.5: Molecular docking module.

Docks candidates against the target protein using AutoDock Vina.
Runs between M4 (Structure Prediction) and M5 (Binding & Scoring)
so that binding_score values are populated before composite scoring.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from drugdiscovery.types import Candidate, PipelineConfig, TargetProfile
from drugdiscovery.utils.io import write_csv

logger = logging.getLogger(__name__)


def _compute_centroid_from_pdb(
    pdb_path: str,
    residue_numbers: Optional[list[int]] = None,
) -> Optional[tuple[float, float, float]]:
    """Compute docking center from receptor PDB alpha-carbon coordinates.

    If *residue_numbers* is provided, averages CA atoms of those residues only
    (i.e. the binding pocket).  Otherwise averages all CA atoms (protein center).

    Returns None if no CA atoms are found.
    """
    try:
        from Bio.PDB import PDBParser
    except ImportError:
        logger.warning("Biopython required for PDB centroid calculation")
        return None

    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("receptor", pdb_path)
    except Exception as exc:
        logger.warning("Failed to parse receptor PDB %s: %s", pdb_path, exc)
        return None

    xs, ys, zs = [], [], []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":
                    continue  # skip hetero atoms
                if residue_numbers and residue.id[1] not in residue_numbers:
                    continue
                if "CA" in residue:
                    coord = residue["CA"].get_vector()
                    xs.append(coord[0])
                    ys.append(coord[1])
                    zs.append(coord[2])
        break  # first model only

    if not xs:
        return None

    n = len(xs)
    return (
        round(sum(xs) / n, 3),
        round(sum(ys) / n, 3),
        round(sum(zs) / n, 3),
    )


def dock_candidates(
    cfg: PipelineConfig,
    candidates: list[Candidate],
    target_profile: TargetProfile,
    run_dir: Path,
) -> list[Candidate]:
    """Dock all candidates against the target protein receptor.

    Populates ``candidate.binding_score`` (kcal/mol, negative = good).
    If the receptor PDB is unavailable or docking tools are missing,
    candidates are returned unchanged (graceful skip).
    """
    receptor_pdb = target_profile.structure_pdb_path
    if not receptor_pdb or not Path(receptor_pdb).exists():
        logger.warning(
            "[M4.5] No receptor PDB available (structure_pdb_path=%s); "
            "skipping docking -- binding_score stays at 0.0",
            receptor_pdb,
        )
        return candidates

    # --- Determine docking center ---
    center: Optional[tuple[float, float, float]] = None

    # Try binding pocket centroid first
    for pocket in target_profile.binding_pockets:
        if pocket.centroid:
            center = pocket.centroid
            logger.info("[M4.5] Using pocket centroid: %s", center)
            break

    # Fallback: compute from PDB
    if center is None:
        pocket_residues = None
        for pocket in target_profile.binding_pockets:
            if pocket.residue_numbers:
                pocket_residues = pocket.residue_numbers
                break
        center = _compute_centroid_from_pdb(receptor_pdb, pocket_residues)
        if center is None:
            logger.warning("[M4.5] Could not determine docking center; skipping docking")
            return candidates
        logger.info("[M4.5] Computed docking center from PDB: %s", center)

    # --- Docking parameters ---
    vina_cfg = cfg.vina_settings
    exhaustiveness = vina_cfg.get("exhaustiveness", 32)
    n_poses = vina_cfg.get("n_poses", 10)
    box_size = vina_cfg.get("box_size", 30.0)

    # --- Import vina tools ---
    try:
        from drugdiscovery.tools.vina import dock_smiles, dock_peptide_pdb
    except ImportError:
        logger.warning("[M4.5] Vina tools unavailable; skipping docking")
        return candidates

    # --- Dock each candidate ---
    output_dir = run_dir / "docking"
    output_dir.mkdir(parents=True, exist_ok=True)

    docked_count = 0
    failed_count = 0
    scores: list[float] = []
    result_rows: list[dict] = []

    for i, cand in enumerate(candidates):
        logger.debug(
            "[M4.5] Docking candidate %d/%d: %s (%s)",
            i + 1, len(candidates), cand.candidate_id, cand.modality,
        )
        score: Optional[float] = None

        if cand.modality == "small_molecule" and cand.smiles:
            score = dock_smiles(
                smiles=cand.smiles,
                receptor_pdb=receptor_pdb,
                center=center,
                box_size=box_size,
                exhaustiveness=exhaustiveness,
                n_poses=n_poses,
            )

        elif cand.modality == "peptide":
            peptide_pdb = (
                cand.metadata.get("conformer_pdb")
                or cand.metadata.get("structure_pdb")
                or cand.metadata.get("complex_pdb")
            )
            if peptide_pdb and Path(peptide_pdb).exists():
                score = dock_peptide_pdb(
                    peptide_pdb=peptide_pdb,
                    receptor_pdb=receptor_pdb,
                    center=center,
                    box_size=box_size,
                    exhaustiveness=exhaustiveness,
                )
            else:
                logger.debug(
                    "[M4.5] No PDB structure for peptide %s; skipping dock",
                    cand.candidate_id,
                )

        if score is not None:
            cand.binding_score = score
            docked_count += 1
            scores.append(score)
        else:
            cand.binding_score = 0.0
            failed_count += 1
            logger.debug("[M4.5] Docking failed for %s; binding_score=0.0", cand.candidate_id)

        result_rows.append({
            "candidate_id": cand.candidate_id,
            "modality": cand.modality,
            "binding_score": round(cand.binding_score, 4),
            "docked": score is not None,
        })

    # --- Save results ---
    if result_rows:
        write_csv(output_dir / "docking_results.csv", result_rows)

    # --- Summary ---
    if scores:
        logger.info(
            "[M4.5] Docking complete: %d/%d docked successfully, %d failed. "
            "Score range: [%.2f, %.2f] kcal/mol",
            docked_count, len(candidates), failed_count,
            min(scores), max(scores),
        )
    else:
        logger.warning(
            "[M4.5] No candidates were successfully docked (%d attempted)",
            len(candidates),
        )

    return candidates
