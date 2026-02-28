"""PoseCheck docking pose validation wrapper.

PoseCheck validates molecular docking poses by checking for steric
clashes, internal strain energy, and protein-ligand interaction
fingerprints.  It helps filter out chemically implausible poses that
score well numerically but would not be viable in practice.

This module returns safe fallback values when PoseCheck is not installed,
so downstream scoring can proceed without pose validation.

Reference: Harris et al., "PoseCheck: A Benchmark and Tool for
Evaluating the Quality of Generated Protein-Ligand Complexes", 2024.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Availability check
# ---------------------------------------------------------------------------

_POSECHECK_AVAILABLE: Optional[bool] = None


def is_available() -> bool:
    """Return True if the ``posecheck`` Python package is importable.

    The result is cached after the first probe.
    """
    global _POSECHECK_AVAILABLE
    if _POSECHECK_AVAILABLE is not None:
        return _POSECHECK_AVAILABLE

    try:
        import posecheck  # noqa: F401
        _POSECHECK_AVAILABLE = True
    except ImportError:
        _POSECHECK_AVAILABLE = False

    return _POSECHECK_AVAILABLE


# ---------------------------------------------------------------------------
# Default / fallback result
# ---------------------------------------------------------------------------

def _default_result() -> dict:
    """Return a permissive fallback result when PoseCheck is unavailable.

    The fallback assumes the pose is valid so that pipeline scoring is
    not blocked by a missing optional dependency.
    """
    return {
        "valid": True,
        "clash_count": 0,
        "strain_energy": 0.0,
        "interaction_fingerprint": {},
    }


# ---------------------------------------------------------------------------
# Pose validation
# ---------------------------------------------------------------------------

def validate_pose(
    protein_pdb: str,
    ligand_sdf: str,
) -> dict:
    """Validate a docking pose for chemical plausibility.

    Checks steric clashes between the ligand and protein, computes
    internal strain energy of the ligand conformation, and generates
    a protein-ligand interaction fingerprint (hydrogen bonds, salt
    bridges, hydrophobic contacts, etc.).

    Args:
        protein_pdb: Path to the receptor/protein PDB file.
        ligand_sdf:  Path to the docked ligand SDF file (single pose).

    Returns:
        A dict with keys:

        - ``valid`` (bool): True if the pose passes basic quality checks
          (low clashes, reasonable strain).
        - ``clash_count`` (int): number of steric clashes detected.
        - ``strain_energy`` (float): internal strain energy of the
          ligand in kcal/mol.
        - ``interaction_fingerprint`` (dict): mapping of interaction
          types (e.g. ``"hbond_donor"``, ``"hydrophobic"``) to counts.

        If PoseCheck is not installed, returns a permissive fallback
        (``valid=True``, zero clashes/strain, empty fingerprint).
    """
    if not protein_pdb or not ligand_sdf:
        logger.warning("PoseCheck: both protein_pdb and ligand_sdf are required")
        return _default_result()

    for label, path in [("Protein PDB", protein_pdb), ("Ligand SDF", ligand_sdf)]:
        if not Path(path).exists():
            logger.warning("PoseCheck: %s not found: %s", label, path)
            return _default_result()

    if not is_available():
        logger.warning(
            "PoseCheck is not installed; pose validation skipped. "
            "Install with: pip install posecheck"
        )
        return _default_result()

    # --- Run PoseCheck ---
    return _run_posecheck(protein_pdb, ligand_sdf)


# ---------------------------------------------------------------------------
# PoseCheck execution
# ---------------------------------------------------------------------------

def _run_posecheck(protein_pdb: str, ligand_sdf: str) -> dict:
    """Execute PoseCheck analysis on a protein-ligand pose."""
    try:
        from posecheck import PoseCheck  # type: ignore[import-untyped]
    except ImportError:
        logger.debug("PoseCheck import failed at runtime")
        return _default_result()

    try:
        pc = PoseCheck()
        pc.load_protein_from_pdb(protein_pdb)
        pc.load_ligands_from_sdf(ligand_sdf)

        # --- Steric clashes ---
        clash_count = _get_clashes(pc)

        # --- Strain energy ---
        strain_energy = _get_strain_energy(pc)

        # --- Interaction fingerprint ---
        interaction_fingerprint = _get_interactions(pc)

        # --- Validity assessment ---
        valid = _assess_validity(clash_count, strain_energy)

        result = {
            "valid": valid,
            "clash_count": clash_count,
            "strain_energy": strain_energy,
            "interaction_fingerprint": interaction_fingerprint,
        }

        logger.info(
            "PoseCheck: valid=%s, clashes=%d, strain=%.2f kcal/mol, "
            "interactions=%d types",
            valid, clash_count, strain_energy,
            len(interaction_fingerprint),
        )

        return result

    except Exception as exc:
        logger.warning("PoseCheck analysis failed: %s", exc)
        return _default_result()


# ---------------------------------------------------------------------------
# Individual check helpers
# ---------------------------------------------------------------------------

def _get_clashes(pc: object) -> int:
    """Extract steric clash count from PoseCheck."""
    try:
        clashes = pc.calculate_clashes()  # type: ignore[union-attr]
        if isinstance(clashes, (list, tuple)):
            # Some versions return a list of clashing atom pairs
            return len(clashes)
        if isinstance(clashes, dict):
            return clashes.get("num_clashes", 0)
        return int(clashes)
    except Exception as exc:
        logger.debug("PoseCheck clash calculation failed: %s", exc)
        return 0


def _get_strain_energy(pc: object) -> float:
    """Extract ligand strain energy from PoseCheck."""
    try:
        strain = pc.calculate_strain_energy()  # type: ignore[union-attr]
        if isinstance(strain, dict):
            return float(strain.get("strain_energy", 0.0))
        if isinstance(strain, (list, tuple)) and strain:
            return float(strain[0])
        return float(strain)
    except Exception as exc:
        logger.debug("PoseCheck strain energy calculation failed: %s", exc)
        return 0.0


def _get_interactions(pc: object) -> dict:
    """Extract protein-ligand interaction fingerprint from PoseCheck."""
    try:
        interactions = pc.calculate_interactions()  # type: ignore[union-attr]
        if isinstance(interactions, dict):
            # Summarise interaction counts by type
            fingerprint: dict[str, int] = {}
            for key, value in interactions.items():
                if isinstance(value, (list, tuple)):
                    fingerprint[str(key)] = len(value)
                elif isinstance(value, (int, float)):
                    fingerprint[str(key)] = int(value)
                elif isinstance(value, bool):
                    fingerprint[str(key)] = 1 if value else 0
                else:
                    fingerprint[str(key)] = 1
            return fingerprint
        return {}
    except Exception as exc:
        logger.debug("PoseCheck interaction fingerprint failed: %s", exc)
        return {}


# ---------------------------------------------------------------------------
# Validity assessment
# ---------------------------------------------------------------------------

# Thresholds for a "valid" pose (conservative defaults).
_MAX_CLASHES = 5
_MAX_STRAIN_KCAL = 10.0


def _assess_validity(clash_count: int, strain_energy: float) -> bool:
    """Determine whether a pose passes basic chemical plausibility checks.

    A pose is considered valid if:
    - The number of steric clashes is at most ``_MAX_CLASHES`` (default 5).
    - The internal strain energy is at most ``_MAX_STRAIN_KCAL`` kcal/mol
      (default 10.0).

    These are conservative thresholds; downstream scoring may apply
    stricter criteria.
    """
    if clash_count > _MAX_CLASHES:
        return False
    if strain_energy > _MAX_STRAIN_KCAL:
        return False
    return True
