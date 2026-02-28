"""Boltz-1 structure and affinity prediction wrapper.

Boltz-1 is an open-source model for predicting protein-ligand complex
structures, comparable to AlphaFold3 for biomolecular complexes.  It
predicts 3D coordinates and can estimate binding affinity from the
predicted complex.

This module is GPU-intensive; CUDA availability is checked before
attempting predictions.  Falls back gracefully when Boltz is not
installed.

Reference: Wohlwend et al., "Boltz-1: Democratizing Biomolecular
Interaction Modeling", 2024.  MIT License.
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

_BOLTZ_AVAILABLE: Optional[bool] = None


def is_available() -> bool:
    """Return True if Boltz-1 is usable (Python package or CLI).

    Checks for the ``boltz`` Python package first, then falls back to
    probing for a ``boltz`` CLI binary on PATH.  The result is cached.
    """
    global _BOLTZ_AVAILABLE
    if _BOLTZ_AVAILABLE is not None:
        return _BOLTZ_AVAILABLE

    # 1. Try Python import
    try:
        import boltz  # noqa: F401
        _BOLTZ_AVAILABLE = True
        return True
    except ImportError:
        pass

    # 2. Try CLI binary
    if shutil.which("boltz") is not None:
        _BOLTZ_AVAILABLE = True
        return True

    _BOLTZ_AVAILABLE = False
    return False


def _cuda_available() -> bool:
    """Check whether a CUDA-capable GPU is accessible via PyTorch."""
    try:
        import torch
        return torch.cuda.is_available()
    except ImportError:
        return False


# ---------------------------------------------------------------------------
# Complex prediction
# ---------------------------------------------------------------------------

def predict_complex(
    protein_sequence: str,
    ligand_smiles: str,
    output_dir: str,
) -> Optional[dict]:
    """Predict a protein-ligand complex structure using Boltz-1.

    Args:
        protein_sequence: Amino acid sequence of the target protein.
        ligand_smiles:    SMILES string of the small-molecule ligand.
        output_dir:       Directory to write output files (PDB, scores).

    Returns:
        A dict with keys:

        - ``pdb_path`` (str): path to the predicted complex PDB file.
        - ``confidence`` (float): model confidence score (0-1).
        - ``predicted_affinity`` (float): estimated binding affinity
          in kcal/mol, or ``NaN`` if not available.

        Returns ``None`` if Boltz is not installed, no GPU is available,
        or prediction fails.
    """
    if not protein_sequence or not ligand_smiles:
        logger.warning("Boltz: protein_sequence and ligand_smiles are both required")
        return None

    if not is_available():
        logger.warning(
            "Boltz-1 is not installed; structure prediction unavailable. "
            "Install with: pip install boltz"
        )
        return None

    if not _cuda_available():
        logger.warning(
            "Boltz-1 requires a CUDA-capable GPU but none was detected. "
            "Skipping complex prediction."
        )
        return None

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # --- Try Python API first ---
    result = _predict_api(protein_sequence, ligand_smiles, out)
    if result is not None:
        return result

    # --- Fallback: CLI ---
    result = _predict_cli(protein_sequence, ligand_smiles, out)
    if result is not None:
        return result

    logger.warning("Boltz-1 prediction failed for the given inputs")
    return None


# ---------------------------------------------------------------------------
# Python API backend
# ---------------------------------------------------------------------------

def _predict_api(
    protein_sequence: str,
    ligand_smiles: str,
    output_dir: Path,
) -> Optional[dict]:
    """Attempt Boltz-1 prediction via the Python API."""
    try:
        from boltz import Boltz1  # type: ignore[import-untyped]
    except ImportError:
        return None

    try:
        model = Boltz1.load()

        prediction = model.predict(
            sequences={
                "protein": protein_sequence,
                "ligand": ligand_smiles,
            },
            output_dir=str(output_dir),
        )

        # Extract results from prediction object
        pdb_path = ""
        confidence = 0.0
        predicted_affinity = float("nan")

        if hasattr(prediction, "output_path"):
            pdb_path = str(prediction.output_path)
        elif hasattr(prediction, "pdb_path"):
            pdb_path = str(prediction.pdb_path)
        else:
            # Search output directory for generated PDB files
            pdb_files = sorted(output_dir.glob("*.pdb"))
            if pdb_files:
                pdb_path = str(pdb_files[0])

        if hasattr(prediction, "confidence"):
            confidence = float(prediction.confidence)
        elif hasattr(prediction, "plddt"):
            confidence = float(prediction.plddt) / 100.0

        if hasattr(prediction, "affinity"):
            predicted_affinity = float(prediction.affinity)
        elif hasattr(prediction, "predicted_affinity"):
            predicted_affinity = float(prediction.predicted_affinity)

        logger.info(
            "Boltz-1 API: confidence=%.3f, affinity=%.3f, pdb=%s",
            confidence, predicted_affinity, pdb_path,
        )

        return {
            "pdb_path": pdb_path,
            "confidence": confidence,
            "predicted_affinity": predicted_affinity,
        }

    except Exception as exc:
        logger.debug("Boltz-1 Python API failed: %s", exc)
        return None


# ---------------------------------------------------------------------------
# CLI backend
# ---------------------------------------------------------------------------

def _predict_cli(
    protein_sequence: str,
    ligand_smiles: str,
    output_dir: Path,
) -> Optional[dict]:
    """Attempt Boltz-1 prediction via the CLI binary."""
    if shutil.which("boltz") is None:
        return None

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)

        # Write input YAML file (Boltz-1 input format)
        input_yaml = tmp / "input.yaml"
        input_content = (
            "sequences:\n"
            "  - protein:\n"
            f"      id: A\n"
            f"      sequence: {protein_sequence}\n"
            "  - ligand:\n"
            f"      id: B\n"
            f"      smiles: {ligand_smiles}\n"
        )
        input_yaml.write_text(input_content, encoding="utf-8")

        cmd = [
            "boltz", "predict",
            str(input_yaml),
            "--out_dir", str(output_dir),
            "--use_msa_server",
        ]

        try:
            proc = subprocess.run(
                cmd, capture_output=True, text=True, timeout=3600,
            )
            if proc.returncode != 0:
                logger.warning("Boltz CLI failed (exit %d): %s", proc.returncode, proc.stderr[:300])
                return None

            return _parse_boltz_output(output_dir)

        except FileNotFoundError:
            logger.warning("boltz binary not found in PATH")
            return None
        except subprocess.TimeoutExpired:
            logger.warning("Boltz CLI timed out (60 min)")
            return None


# ---------------------------------------------------------------------------
# Output parsing
# ---------------------------------------------------------------------------

def _parse_boltz_output(output_dir: Path) -> Optional[dict]:
    """Parse Boltz-1 output directory for the predicted complex.

    Looks for PDB output files and associated confidence/score files.
    """
    import json

    pdb_path = ""
    confidence = 0.0
    predicted_affinity = float("nan")

    # Find PDB output
    pdb_files = sorted(output_dir.rglob("*.pdb"))
    if pdb_files:
        pdb_path = str(pdb_files[0])

    if not pdb_path:
        cif_files = sorted(output_dir.rglob("*.cif"))
        if cif_files:
            pdb_path = str(cif_files[0])

    if not pdb_path:
        logger.warning("No structure output found in Boltz output directory")
        return None

    # Parse confidence/scores from JSON files
    for json_path in output_dir.rglob("*.json"):
        try:
            data = json.loads(json_path.read_text(encoding="utf-8"))

            if "confidence" in data:
                confidence = float(data["confidence"])
            elif "plddt" in data:
                plddt_vals = data["plddt"]
                if isinstance(plddt_vals, list) and plddt_vals:
                    confidence = sum(plddt_vals) / len(plddt_vals) / 100.0
                elif isinstance(plddt_vals, (int, float)):
                    confidence = float(plddt_vals) / 100.0

            if "affinity" in data:
                predicted_affinity = float(data["affinity"])
            elif "predicted_affinity" in data:
                predicted_affinity = float(data["predicted_affinity"])
            elif "binding_energy" in data:
                predicted_affinity = float(data["binding_energy"])

        except Exception as exc:
            logger.debug("Failed to parse Boltz score file %s: %s", json_path, exc)
            continue

    logger.info(
        "Boltz-1 CLI: confidence=%.3f, affinity=%.3f, pdb=%s",
        confidence, predicted_affinity, pdb_path,
    )

    return {
        "pdb_path": pdb_path,
        "confidence": confidence,
        "predicted_affinity": predicted_affinity,
    }
