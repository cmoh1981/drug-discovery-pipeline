"""DiffDock ML blind docking wrapper.

DiffDock is a diffusion-based generative model for blind molecular
docking that does not require a pre-defined binding pocket.  This module
wraps both the Python API (if installed) and the CLI fallback, and
degrades gracefully when DiffDock is unavailable.

Reference: Corso et al., "DiffDock: Diffusion Steps, Twists, and Turns
for Molecular Docking", ICLR 2023.  MIT License.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Availability check
# ---------------------------------------------------------------------------

_DIFFDOCK_AVAILABLE: Optional[bool] = None


def is_available() -> bool:
    """Return True if DiffDock is usable (Python package or CLI).

    The result is cached after the first call to avoid repeated
    import / PATH probes on every invocation.
    """
    global _DIFFDOCK_AVAILABLE
    if _DIFFDOCK_AVAILABLE is not None:
        return _DIFFDOCK_AVAILABLE

    # 1. Try Python import
    try:
        import diffdock  # noqa: F401
        _DIFFDOCK_AVAILABLE = True
        return True
    except ImportError:
        pass

    # 2. Try CLI binary
    if shutil.which("diffdock") is not None:
        _DIFFDOCK_AVAILABLE = True
        return True

    _DIFFDOCK_AVAILABLE = False
    return False


# ---------------------------------------------------------------------------
# Small-molecule docking
# ---------------------------------------------------------------------------

def dock_smiles(
    smiles: str,
    receptor_pdb: str,
    n_poses: int = 10,
) -> Optional[float]:
    """Dock a small molecule (SMILES) against a receptor PDB using DiffDock.

    DiffDock performs *blind* docking -- no binding-site coordinates are
    required.  The best (lowest) predicted binding score in kcal/mol is
    returned.

    Args:
        smiles:       SMILES string of the ligand.
        receptor_pdb: Path to the receptor PDB file.
        n_poses:      Number of poses to sample (default 10).

    Returns:
        Best binding affinity in kcal/mol (negative = favorable),
        or ``None`` if DiffDock is not installed or docking fails.
    """
    if not smiles or not receptor_pdb:
        return None

    if not Path(receptor_pdb).exists():
        logger.warning("Receptor PDB not found: %s", receptor_pdb)
        return None

    # --- Try Python API first ---
    result = _dock_smiles_api(smiles, receptor_pdb, n_poses)
    if result is not None:
        return result

    # --- Fallback: CLI ---
    result = _dock_smiles_cli(smiles, receptor_pdb, n_poses)
    if result is not None:
        return result

    if not is_available():
        logger.warning(
            "DiffDock is not installed; blind docking unavailable. "
            "Install with: pip install diffdock  (or place the diffdock binary on PATH)"
        )
    return None


def dock_peptide_pdb(
    peptide_pdb: str,
    receptor_pdb: str,
    n_poses: int = 10,
) -> Optional[float]:
    """Dock a peptide PDB against a receptor PDB using DiffDock.

    Args:
        peptide_pdb:  Path to the peptide/ligand PDB file.
        receptor_pdb: Path to the receptor PDB file.
        n_poses:      Number of poses to sample.

    Returns:
        Best binding affinity in kcal/mol, or ``None`` on failure.
    """
    if not peptide_pdb or not receptor_pdb:
        return None

    for label, path in [("Peptide PDB", peptide_pdb), ("Receptor PDB", receptor_pdb)]:
        if not Path(path).exists():
            logger.warning("%s not found: %s", label, path)
            return None

    # --- Try Python API ---
    result = _dock_peptide_api(peptide_pdb, receptor_pdb, n_poses)
    if result is not None:
        return result

    # --- Fallback: CLI ---
    result = _dock_peptide_cli(peptide_pdb, receptor_pdb, n_poses)
    if result is not None:
        return result

    if not is_available():
        logger.warning(
            "DiffDock is not installed; peptide blind docking unavailable."
        )
    return None


# ---------------------------------------------------------------------------
# Python API backend
# ---------------------------------------------------------------------------

def _dock_smiles_api(
    smiles: str,
    receptor_pdb: str,
    n_poses: int,
) -> Optional[float]:
    """Attempt DiffDock docking via the Python API."""
    try:
        from diffdock import DiffDockModel  # type: ignore[import-untyped]
    except ImportError:
        return None

    try:
        model = DiffDockModel()
        results = model.dock(
            protein_path=receptor_pdb,
            ligand=smiles,
            num_poses=n_poses,
        )
        if results and hasattr(results[0], "score"):
            scores = [r.score for r in results if r.score is not None]
            if scores:
                best = min(scores)
                logger.info(
                    "DiffDock API: best score %.3f kcal/mol from %d poses",
                    best, len(scores),
                )
                return float(best)
    except Exception as exc:
        logger.debug("DiffDock Python API failed: %s", exc)

    return None


def _dock_peptide_api(
    peptide_pdb: str,
    receptor_pdb: str,
    n_poses: int,
) -> Optional[float]:
    """Attempt DiffDock peptide docking via the Python API."""
    try:
        from diffdock import DiffDockModel  # type: ignore[import-untyped]
    except ImportError:
        return None

    try:
        model = DiffDockModel()
        results = model.dock(
            protein_path=receptor_pdb,
            ligand=peptide_pdb,
            num_poses=n_poses,
        )
        if results and hasattr(results[0], "score"):
            scores = [r.score for r in results if r.score is not None]
            if scores:
                best = min(scores)
                logger.info(
                    "DiffDock API (peptide): best score %.3f kcal/mol from %d poses",
                    best, len(scores),
                )
                return float(best)
    except Exception as exc:
        logger.debug("DiffDock Python API (peptide) failed: %s", exc)

    return None


# ---------------------------------------------------------------------------
# CLI backend
# ---------------------------------------------------------------------------

def _dock_smiles_cli(
    smiles: str,
    receptor_pdb: str,
    n_poses: int,
) -> Optional[float]:
    """Attempt DiffDock docking via the CLI binary."""
    if shutil.which("diffdock") is None:
        return None

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)

        # Write SMILES to a ligand file
        ligand_smi = tmp / "ligand.smi"
        ligand_smi.write_text(smiles, encoding="utf-8")

        output_dir = tmp / "output"
        output_dir.mkdir()

        cmd = [
            "diffdock",
            "--protein_path", str(receptor_pdb),
            "--ligand", str(ligand_smi),
            "--out_dir", str(output_dir),
            "--samples_per_complex", str(n_poses),
        ]

        try:
            proc = subprocess.run(
                cmd, capture_output=True, text=True, timeout=1200,
            )
            if proc.returncode != 0:
                logger.warning("DiffDock CLI failed: %s", proc.stderr[:300])
                return None

            return _parse_diffdock_output(output_dir)

        except FileNotFoundError:
            logger.warning("diffdock binary not found in PATH")
            return None
        except subprocess.TimeoutExpired:
            logger.warning("DiffDock CLI timed out (20 min)")
            return None


def _dock_peptide_cli(
    peptide_pdb: str,
    receptor_pdb: str,
    n_poses: int,
) -> Optional[float]:
    """Attempt DiffDock peptide docking via the CLI binary."""
    if shutil.which("diffdock") is None:
        return None

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        output_dir = tmp / "output"
        output_dir.mkdir()

        cmd = [
            "diffdock",
            "--protein_path", str(receptor_pdb),
            "--ligand", str(peptide_pdb),
            "--out_dir", str(output_dir),
            "--samples_per_complex", str(n_poses),
        ]

        try:
            proc = subprocess.run(
                cmd, capture_output=True, text=True, timeout=1200,
            )
            if proc.returncode != 0:
                logger.warning("DiffDock CLI (peptide) failed: %s", proc.stderr[:300])
                return None

            return _parse_diffdock_output(output_dir)

        except FileNotFoundError:
            logger.warning("diffdock binary not found in PATH")
            return None
        except subprocess.TimeoutExpired:
            logger.warning("DiffDock CLI (peptide) timed out (20 min)")
            return None


# ---------------------------------------------------------------------------
# Output parsing
# ---------------------------------------------------------------------------

def _parse_diffdock_output(output_dir: Path) -> Optional[float]:
    """Parse DiffDock output directory for the best binding score.

    DiffDock typically writes ranked PDB/SDF files whose filenames contain
    the confidence score (e.g., ``rank1_confidence-1.42.sdf``).  It may
    also write a ``scores.csv`` or embed scores as remarks in the SDF.
    """
    import re

    best_score: Optional[float] = None

    # Strategy 1: scores CSV file
    for csv_path in output_dir.rglob("*.csv"):
        try:
            for line in csv_path.read_text(encoding="utf-8").splitlines()[1:]:
                parts = line.strip().split(",")
                for part in parts:
                    try:
                        val = float(part)
                        if best_score is None or val < best_score:
                            best_score = val
                    except ValueError:
                        continue
        except Exception:
            continue

    if best_score is not None:
        return best_score

    # Strategy 2: score embedded in filename (rank1_confidence-1.42.sdf)
    pattern = re.compile(r"confidence[_-]?([-+]?\d+\.?\d*)")
    for f in output_dir.rglob("*rank*"):
        match = pattern.search(f.name)
        if match:
            try:
                val = float(match.group(1))
                if best_score is None or val < best_score:
                    best_score = val
            except ValueError:
                continue

    if best_score is not None:
        return best_score

    # Strategy 3: parse SDF remark lines
    for sdf_path in output_dir.rglob("*.sdf"):
        try:
            for line in sdf_path.read_text(encoding="utf-8").splitlines():
                if "minimizedAffinity" in line or "score" in line.lower():
                    nums = re.findall(r"[-+]?\d+\.?\d*", line)
                    for n in nums:
                        val = float(n)
                        if best_score is None or val < best_score:
                            best_score = val
        except Exception:
            continue

    return best_score
