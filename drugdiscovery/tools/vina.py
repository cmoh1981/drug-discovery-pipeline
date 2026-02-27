"""AutoDock Vina docking wrapper."""

from __future__ import annotations

import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


def dock_smiles(
    smiles: str,
    receptor_pdb: str,
    center: tuple[float, float, float],
    box_size: float = 30.0,
    exhaustiveness: int = 32,
    n_poses: int = 10,
) -> Optional[float]:
    """Dock a small molecule (SMILES) against a receptor PDB.

    Returns best binding affinity in kcal/mol, or None if docking fails.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        logger.warning("RDKit required for ligand preparation")
        return None

    # Prepare ligand
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if result != 0:
        return None
    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        ligand_pdb = tmp / "ligand.pdb"
        Chem.MolToPDBFile(mol, str(ligand_pdb))

        # Try meeko for PDBQT conversion
        ligand_pdbqt = tmp / "ligand.pdbqt"
        receptor_pdbqt = tmp / "receptor.pdbqt"

        if not _prepare_pdbqt(ligand_pdb, ligand_pdbqt, is_ligand=True):
            return None
        if not _prepare_pdbqt(Path(receptor_pdb), receptor_pdbqt, is_ligand=False):
            return None

        # Run Vina
        output_pdbqt = tmp / "output.pdbqt"
        try:
            cmd = [
                "vina",
                "--receptor", str(receptor_pdbqt),
                "--ligand", str(ligand_pdbqt),
                "--center_x", str(center[0]),
                "--center_y", str(center[1]),
                "--center_z", str(center[2]),
                "--size_x", str(box_size),
                "--size_y", str(box_size),
                "--size_z", str(box_size),
                "--exhaustiveness", str(exhaustiveness),
                "--num_modes", str(n_poses),
                "--out", str(output_pdbqt),
            ]
            proc = subprocess.run(
                cmd, capture_output=True, text=True, timeout=600
            )
            if proc.returncode == 0:
                return _parse_vina_output(proc.stdout)
            else:
                logger.warning("Vina failed: %s", proc.stderr[:200])
                return None
        except FileNotFoundError:
            logger.warning("Vina executable not found in PATH")
            return None
        except subprocess.TimeoutExpired:
            logger.warning("Vina timed out")
            return None


def dock_peptide_pdb(
    peptide_pdb: str,
    receptor_pdb: str,
    center: tuple[float, float, float],
    box_size: float = 30.0,
    exhaustiveness: int = 32,
) -> Optional[float]:
    """Dock a peptide PDB against a receptor. Returns best affinity."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        ligand_pdbqt = tmp / "peptide.pdbqt"
        receptor_pdbqt = tmp / "receptor.pdbqt"

        if not _prepare_pdbqt(Path(peptide_pdb), ligand_pdbqt, is_ligand=True):
            return None
        if not _prepare_pdbqt(Path(receptor_pdb), receptor_pdbqt, is_ligand=False):
            return None

        output_pdbqt = tmp / "output.pdbqt"
        try:
            cmd = [
                "vina",
                "--receptor", str(receptor_pdbqt),
                "--ligand", str(ligand_pdbqt),
                "--center_x", str(center[0]),
                "--center_y", str(center[1]),
                "--center_z", str(center[2]),
                "--size_x", str(box_size),
                "--size_y", str(box_size),
                "--size_z", str(box_size),
                "--exhaustiveness", str(exhaustiveness),
                "--out", str(output_pdbqt),
            ]
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            if proc.returncode == 0:
                return _parse_vina_output(proc.stdout)
            return None
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return None


def _prepare_pdbqt(input_path: Path, output_path: Path, is_ligand: bool) -> bool:
    """Convert PDB to PDBQT using meeko (ligand) or obabel (receptor)."""
    if is_ligand:
        try:
            from meeko import MoleculePreparation, PDBQTWriterLegacy
            from rdkit import Chem

            mol = Chem.MolFromPDBFile(str(input_path), removeHs=False)
            if mol is None:
                return False
            prep = MoleculePreparation()
            mol_setup = prep.prepare(mol)[0]
            pdbqt_str, _, _ = PDBQTWriterLegacy.write_string(mol_setup)
            output_path.write_text(pdbqt_str)
            return True
        except ImportError:
            pass
        except Exception as exc:
            logger.debug("Meeko preparation failed: %s", exc)

    # Fallback: try obabel
    try:
        cmd = ["obabel", str(input_path), "-O", str(output_path)]
        if not is_ligand:
            cmd.append("-xr")
        proc = subprocess.run(
            cmd,
            capture_output=True, text=True, timeout=60,
        )
        return proc.returncode == 0 and output_path.exists()
    except FileNotFoundError:
        logger.debug("obabel not found")
        return False


def _parse_vina_output(stdout: str) -> Optional[float]:
    """Parse Vina stdout for best binding affinity."""
    for line in stdout.splitlines():
        line = line.strip()
        if line.startswith("1"):
            parts = line.split()
            if len(parts) >= 2:
                try:
                    return float(parts[1])
                except ValueError:
                    pass
    return None
