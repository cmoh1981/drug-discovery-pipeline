"""GNINA CNN-enhanced molecular docking wrapper.

GNINA is a fork of AutoDock Vina that adds convolutional neural network
(CNN) scoring to traditional docking.  It uses the same CLI interface as
Vina but produces CNN-refined affinity predictions.

This module checks for the ``gnina`` binary on PATH and falls back
gracefully when it is not installed.

Reference: McNutt et al., "GNINA 1.0: molecular docking with deep
learning", J. Cheminform. 2021.
"""

from __future__ import annotations

import logging
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Availability check
# ---------------------------------------------------------------------------

_GNINA_AVAILABLE: Optional[bool] = None


def is_available() -> bool:
    """Return True if the ``gnina`` binary is found on PATH.

    The result is cached after the first probe.
    """
    global _GNINA_AVAILABLE
    if _GNINA_AVAILABLE is not None:
        return _GNINA_AVAILABLE
    _GNINA_AVAILABLE = shutil.which("gnina") is not None
    return _GNINA_AVAILABLE


# ---------------------------------------------------------------------------
# Small-molecule docking
# ---------------------------------------------------------------------------

def dock_smiles(
    smiles: str,
    receptor_pdb: str,
    center: tuple[float, float, float],
    box_size: float = 30.0,
    exhaustiveness: int = 32,
    n_poses: int = 10,
) -> Optional[float]:
    """Dock a small molecule (SMILES) against a receptor PDB using GNINA.

    The ligand is converted from SMILES to SDF via RDKit, then docked
    with GNINA's CNN-rescoring pipeline.

    Args:
        smiles:          SMILES string of the ligand.
        receptor_pdb:    Path to the receptor PDB file.
        center:          (x, y, z) center of the docking box in Angstroms.
        box_size:        Side length of the cubic docking box (default 30 A).
        exhaustiveness:  Search exhaustiveness (default 32).
        n_poses:         Maximum number of poses to generate (default 10).

    Returns:
        Best CNN-rescored affinity in kcal/mol (negative = favorable),
        or ``None`` if GNINA is not installed or docking fails.
    """
    if not smiles or not receptor_pdb:
        return None

    if not Path(receptor_pdb).exists():
        logger.warning("Receptor PDB not found: %s", receptor_pdb)
        return None

    if not is_available():
        logger.warning(
            "GNINA binary not found on PATH; CNN-enhanced docking unavailable. "
            "See https://github.com/gnina/gnina for installation."
        )
        return None

    # --- Prepare ligand SDF with RDKit ---
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        logger.warning("RDKit required for SMILES-to-SDF ligand preparation")
        return None

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.warning("RDKit failed to parse SMILES: %.60s", smiles)
        return None

    mol = Chem.AddHs(mol)
    embed_result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if embed_result != 0:
        logger.warning("RDKit 3D embedding failed for SMILES: %.60s", smiles)
        return None
    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        ligand_sdf = tmp / "ligand.sdf"
        writer = Chem.SDWriter(str(ligand_sdf))
        writer.write(mol)
        writer.close()

        return _run_gnina(
            receptor_pdb=receptor_pdb,
            ligand_path=str(ligand_sdf),
            center=center,
            box_size=box_size,
            exhaustiveness=exhaustiveness,
            n_poses=n_poses,
            output_dir=tmp,
        )


# ---------------------------------------------------------------------------
# Peptide docking
# ---------------------------------------------------------------------------

def dock_peptide_pdb(
    peptide_pdb: str,
    receptor_pdb: str,
    center: tuple[float, float, float],
    box_size: float = 30.0,
) -> Optional[float]:
    """Dock a peptide PDB against a receptor PDB using GNINA.

    Args:
        peptide_pdb:  Path to the peptide/ligand PDB file.
        receptor_pdb: Path to the receptor PDB file.
        center:       (x, y, z) center of the docking box in Angstroms.
        box_size:     Side length of the cubic docking box (default 30 A).

    Returns:
        Best CNN-rescored affinity in kcal/mol, or ``None`` on failure.
    """
    if not peptide_pdb or not receptor_pdb:
        return None

    for label, path in [("Peptide PDB", peptide_pdb), ("Receptor PDB", receptor_pdb)]:
        if not Path(path).exists():
            logger.warning("%s not found: %s", label, path)
            return None

    if not is_available():
        logger.warning("GNINA binary not found on PATH; peptide docking unavailable.")
        return None

    with tempfile.TemporaryDirectory() as tmpdir:
        return _run_gnina(
            receptor_pdb=receptor_pdb,
            ligand_path=peptide_pdb,
            center=center,
            box_size=box_size,
            exhaustiveness=32,
            n_poses=10,
            output_dir=Path(tmpdir),
        )


# ---------------------------------------------------------------------------
# Core GNINA execution
# ---------------------------------------------------------------------------

def _run_gnina(
    receptor_pdb: str,
    ligand_path: str,
    center: tuple[float, float, float],
    box_size: float,
    exhaustiveness: int,
    n_poses: int,
    output_dir: Path,
) -> Optional[float]:
    """Execute the ``gnina`` binary and return the best affinity score.

    GNINA uses the same CLI flags as AutoDock Vina, with additional
    CNN-specific options enabled by default.
    """
    output_sdf = output_dir / "gnina_output.sdf"

    cmd = [
        "gnina",
        "--receptor", str(receptor_pdb),
        "--ligand", str(ligand_path),
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(box_size),
        "--size_y", str(box_size),
        "--size_z", str(box_size),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(n_poses),
        "--out", str(output_sdf),
        "--cnn_scoring", "rescore",
    ]

    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=600,
        )
        if proc.returncode != 0:
            logger.warning("GNINA failed (exit %d): %s", proc.returncode, proc.stderr[:300])
            return None

        score = _parse_gnina_stdout(proc.stdout)
        if score is not None:
            logger.info("GNINA: best affinity %.3f kcal/mol", score)
        return score

    except FileNotFoundError:
        logger.warning("GNINA binary not found in PATH")
        return None
    except subprocess.TimeoutExpired:
        logger.warning("GNINA timed out (10 min)")
        return None


# ---------------------------------------------------------------------------
# Output parsing
# ---------------------------------------------------------------------------

def _parse_gnina_stdout(stdout: str) -> Optional[float]:
    """Parse GNINA stdout for the best binding affinity.

    GNINA output follows the same tabular format as Vina::

        mode |   affinity | dist from best mode
             | (kcal/mol) | rmsd l.b.| rmsd u.b.
        -----+------------+----------+----------
           1       -8.320      0.000      0.000
           2       -7.450      2.134      3.567
           ...

    The first numeric mode line gives the best (most negative) score.
    """
    in_results = False
    for line in stdout.splitlines():
        stripped = line.strip()

        # Detect the separator line before results
        if stripped.startswith("-----"):
            in_results = True
            continue

        if in_results and stripped:
            parts = stripped.split()
            if len(parts) >= 2:
                try:
                    return float(parts[1])
                except ValueError:
                    pass

    # Fallback: scan for any line starting with "1" (first pose)
    for line in stdout.splitlines():
        stripped = line.strip()
        if re.match(r"^\s*1\s+[-]?\d+\.\d+", stripped):
            parts = stripped.split()
            if len(parts) >= 2:
                try:
                    return float(parts[1])
                except ValueError:
                    pass

    return None
