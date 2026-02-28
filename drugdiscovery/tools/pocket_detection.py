"""Binding pocket detection wrapper (P2Rank / fpocket / fallback).

Automatically detects binding pockets from a PDB structure file.
Supports P2Rank (Java CLI) and fpocket (C CLI) with graceful
fallback to centroid-based pocket estimation when neither is installed.

References:
  - P2Rank: Krivak & Hoksza, J. Cheminform. 2018.
  - fpocket: Le Guilloux et al., BMC Bioinformatics, 2009.
"""

from __future__ import annotations

import logging
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

from drugdiscovery.types import BindingPocket

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Availability checks
# ---------------------------------------------------------------------------

_P2RANK_AVAILABLE: Optional[bool] = None
_FPOCKET_AVAILABLE: Optional[bool] = None


def _check_p2rank() -> bool:
    """Return True if the ``prank`` (P2Rank) binary is on PATH."""
    global _P2RANK_AVAILABLE
    if _P2RANK_AVAILABLE is not None:
        return _P2RANK_AVAILABLE
    _P2RANK_AVAILABLE = shutil.which("prank") is not None
    return _P2RANK_AVAILABLE


def _check_fpocket() -> bool:
    """Return True if the ``fpocket`` binary is on PATH."""
    global _FPOCKET_AVAILABLE
    if _FPOCKET_AVAILABLE is not None:
        return _FPOCKET_AVAILABLE
    _FPOCKET_AVAILABLE = shutil.which("fpocket") is not None
    return _FPOCKET_AVAILABLE


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def detect_pockets(pdb_path: str, max_pockets: int = 5) -> list[BindingPocket]:
    """Detect binding pockets from a PDB structure file.

    Tries the following tools in order:
      1. P2Rank (``prank predict``)
      2. fpocket (``fpocket -f``)
      3. Biopython CA centroid fallback (whole protein center)

    Args:
        pdb_path:    Path to the receptor PDB file.
        max_pockets: Maximum number of pockets to return.

    Returns:
        List of BindingPocket instances with centroids populated.
        Returns an empty list only if all methods fail.
    """
    if not pdb_path or not Path(pdb_path).exists():
        logger.warning("Pocket detection: PDB file not found: %s", pdb_path)
        return []

    # Strategy 1: P2Rank
    if _check_p2rank():
        pockets = _detect_p2rank(pdb_path, max_pockets)
        if pockets:
            logger.info("Pocket detection (P2Rank): found %d pockets", len(pockets))
            return pockets

    # Strategy 2: fpocket
    if _check_fpocket():
        pockets = _detect_fpocket(pdb_path, max_pockets)
        if pockets:
            logger.info("Pocket detection (fpocket): found %d pockets", len(pockets))
            return pockets

    # Strategy 3: Biopython centroid fallback
    pockets = _detect_centroid_fallback(pdb_path)
    if pockets:
        logger.info(
            "Pocket detection (centroid fallback): returning whole-protein center. "
            "Install P2Rank or fpocket for better pocket detection."
        )
        return pockets

    logger.warning("Pocket detection: all strategies failed for %s", pdb_path)
    return []


# ---------------------------------------------------------------------------
# P2Rank backend
# ---------------------------------------------------------------------------

def _detect_p2rank(pdb_path: str, max_pockets: int) -> list[BindingPocket]:
    """Run P2Rank (prank predict) and parse results."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        output_dir = tmp / "p2rank_output"

        cmd = [
            "prank", "predict",
            "-f", str(pdb_path),
            "-o", str(output_dir),
        ]

        try:
            proc = subprocess.run(
                cmd, capture_output=True, text=True, timeout=300,
            )
            if proc.returncode != 0:
                logger.warning("P2Rank failed (exit %d): %s", proc.returncode, proc.stderr[:300])
                return []

            return _parse_p2rank_output(output_dir, max_pockets)

        except FileNotFoundError:
            logger.debug("prank binary not found")
            return []
        except subprocess.TimeoutExpired:
            logger.warning("P2Rank timed out (5 min)")
            return []


def _parse_p2rank_output(output_dir: Path, max_pockets: int) -> list[BindingPocket]:
    """Parse P2Rank CSV output for pocket predictions.

    P2Rank writes a ``*_predictions.csv`` with columns:
      name, rank, score, probability, sas_points, surf_atoms,
      center_x, center_y, center_z, residue_ids, surf_atom_ids
    """
    pockets: list[BindingPocket] = []

    csv_files = list(output_dir.rglob("*_predictions.csv"))
    if not csv_files:
        return []

    csv_path = csv_files[0]
    try:
        lines = csv_path.read_text(encoding="utf-8").strip().splitlines()
        if len(lines) < 2:
            return []

        # Parse header to find column indices
        header = [h.strip() for h in lines[0].split(",")]
        col = {name: idx for idx, name in enumerate(header)}

        for line in lines[1 : max_pockets + 1]:
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < len(header):
                continue

            try:
                cx = float(parts[col.get("center_x", 6)])
                cy = float(parts[col.get("center_y", 7)])
                cz = float(parts[col.get("center_z", 8)])
            except (ValueError, IndexError):
                continue

            # Parse residue IDs (format: "A_123 A_124 ...")
            residue_str = parts[col.get("residue_ids", 9)] if "residue_ids" in col else ""
            residue_numbers = _parse_residue_ids(residue_str)

            rank = parts[col.get("rank", 1)] if "rank" in col else str(len(pockets) + 1)
            score = parts[col.get("score", 2)] if "score" in col else "0"

            pockets.append(BindingPocket(
                pocket_id=f"p2rank_pocket_{rank}",
                residue_numbers=residue_numbers,
                description=f"P2Rank pocket (score={score})",
                pocket_type="predicted",
                centroid=(round(cx, 3), round(cy, 3), round(cz, 3)),
            ))

    except Exception as exc:
        logger.debug("Failed to parse P2Rank output: %s", exc)

    return pockets


# ---------------------------------------------------------------------------
# fpocket backend
# ---------------------------------------------------------------------------

def _detect_fpocket(pdb_path: str, max_pockets: int) -> list[BindingPocket]:
    """Run fpocket and parse results."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)

        # fpocket needs the PDB in its working directory
        import shutil as shu
        local_pdb = tmp / Path(pdb_path).name
        shu.copy2(pdb_path, local_pdb)

        cmd = ["fpocket", "-f", str(local_pdb)]

        try:
            proc = subprocess.run(
                cmd, capture_output=True, text=True, timeout=300, cwd=str(tmp),
            )
            if proc.returncode != 0:
                logger.warning("fpocket failed (exit %d): %s", proc.returncode, proc.stderr[:300])
                return []

            return _parse_fpocket_output(tmp, local_pdb.stem, max_pockets)

        except FileNotFoundError:
            logger.debug("fpocket binary not found")
            return []
        except subprocess.TimeoutExpired:
            logger.warning("fpocket timed out (5 min)")
            return []


def _parse_fpocket_output(work_dir: Path, stem: str, max_pockets: int) -> list[BindingPocket]:
    """Parse fpocket output directory.

    fpocket creates ``{stem}_out/`` with pocket PDB files and an info file.
    """
    pockets: list[BindingPocket] = []
    out_dir = work_dir / f"{stem}_out"

    if not out_dir.exists():
        # Try alternative naming
        candidates = list(work_dir.glob("*_out"))
        if candidates:
            out_dir = candidates[0]
        else:
            return []

    # Parse pocket PDB files for centroids
    pocket_pdbs = sorted(out_dir.glob("pockets/pocket*.pdb"))
    if not pocket_pdbs:
        pocket_pdbs = sorted(out_dir.glob("pocket*.pdb"))

    for idx, pocket_pdb in enumerate(pocket_pdbs[:max_pockets], start=1):
        centroid = _compute_pdb_centroid(str(pocket_pdb))
        residues = _extract_residue_numbers_from_pdb(str(pocket_pdb))

        pockets.append(BindingPocket(
            pocket_id=f"fpocket_{idx:02d}",
            residue_numbers=residues,
            description=f"fpocket pocket #{idx}",
            pocket_type="predicted",
            centroid=centroid,
        ))

    return pockets


# ---------------------------------------------------------------------------
# Centroid fallback
# ---------------------------------------------------------------------------

def _detect_centroid_fallback(pdb_path: str) -> list[BindingPocket]:
    """Compute whole-protein CA centroid as a fallback 'pocket'."""
    centroid = _compute_pdb_centroid(pdb_path)
    if centroid is None:
        return []

    return [BindingPocket(
        pocket_id="centroid_fallback",
        residue_numbers=[],
        description="Whole-protein centroid (no pocket detection tool available)",
        pocket_type="centroid",
        centroid=centroid,
    )]


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _compute_pdb_centroid(
    pdb_path: str,
    residue_numbers: Optional[list[int]] = None,
) -> Optional[tuple[float, float, float]]:
    """Compute centroid from CA atoms in a PDB file."""
    try:
        from Bio.PDB import PDBParser
    except ImportError:
        # Pure fallback: parse PDB manually for ATOM lines
        return _compute_centroid_manual(pdb_path, residue_numbers)

    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("pdb", pdb_path)
    except Exception:
        return _compute_centroid_manual(pdb_path, residue_numbers)

    xs, ys, zs = [], [], []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":
                    continue
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


def _compute_centroid_manual(
    pdb_path: str,
    residue_numbers: Optional[list[int]] = None,
) -> Optional[tuple[float, float, float]]:
    """Parse PDB ATOM records manually for CA centroid (no Biopython)."""
    xs, ys, zs = [], [], []
    try:
        with open(pdb_path, encoding="utf-8") as fh:
            for line in fh:
                if not line.startswith(("ATOM  ", "HETATM")):
                    continue
                atom_name = line[12:16].strip()
                if atom_name != "CA":
                    continue
                if residue_numbers:
                    try:
                        resnum = int(line[22:26].strip())
                    except ValueError:
                        continue
                    if resnum not in residue_numbers:
                        continue
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    xs.append(x)
                    ys.append(y)
                    zs.append(z)
                except ValueError:
                    continue
    except Exception:
        return None

    if not xs:
        return None

    n = len(xs)
    return (
        round(sum(xs) / n, 3),
        round(sum(ys) / n, 3),
        round(sum(zs) / n, 3),
    )


def _parse_residue_ids(residue_str: str) -> list[int]:
    """Parse P2Rank residue ID string (e.g. 'A_123 A_124 B_56') to integers."""
    nums: list[int] = []
    for part in residue_str.split():
        match = re.search(r"(\d+)", part)
        if match:
            nums.append(int(match.group(1)))
    return sorted(set(nums))


def _extract_residue_numbers_from_pdb(pdb_path: str) -> list[int]:
    """Extract unique residue numbers from a PDB file."""
    residues: set[int] = set()
    try:
        with open(pdb_path, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith(("ATOM  ", "HETATM")):
                    try:
                        resnum = int(line[22:26].strip())
                        residues.add(resnum)
                    except ValueError:
                        continue
    except Exception:
        pass
    return sorted(residues)
