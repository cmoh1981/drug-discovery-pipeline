"""ColabFold wrapper for protein-peptide complex prediction."""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


def predict_complex(
    target_sequence: str,
    peptide_sequence: str,
    output_dir: str,
    job_name: str = "complex",
    num_models: int = 1,
    num_recycles: int = 3,
) -> Optional[dict]:
    """Predict protein-peptide complex using ColabFold.

    Returns dict with keys: pdb_path, iptm, ptm, plddt
    or None if prediction fails.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Create input FASTA with colon-separated chains
    fasta_path = out / f"{job_name}.fasta"
    fasta_path.write_text(
        f">{job_name}\n{target_sequence}:{peptide_sequence}\n",
        encoding="utf-8",
    )

    try:
        cmd = [
            "colabfold_batch",
            str(fasta_path),
            str(out),
            "--num-models", str(num_models),
            "--num-recycle", str(num_recycles),
            "--model-type", "alphafold2_multimer_v3",
        ]
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=1800
        )

        if proc.returncode != 0:
            logger.warning("ColabFold failed: %s", proc.stderr[:300])
            return None

        # Find output PDB
        pdb_files = list(out.glob(f"{job_name}*rank_001*.pdb"))
        if not pdb_files:
            pdb_files = list(out.glob(f"{job_name}*.pdb"))

        if not pdb_files:
            logger.warning("No PDB output found from ColabFold")
            return None

        pdb_path = pdb_files[0]

        # Parse scores from log or filename
        scores = _parse_colabfold_scores(out, job_name)

        return {
            "pdb_path": str(pdb_path),
            "iptm": scores.get("iptm", 0.0),
            "ptm": scores.get("ptm", 0.0),
            "plddt": scores.get("plddt", 0.0),
        }

    except FileNotFoundError:
        logger.warning("colabfold_batch not found in PATH")
        return None
    except subprocess.TimeoutExpired:
        logger.warning("ColabFold timed out (30 min)")
        return None


def _parse_colabfold_scores(output_dir: Path, job_name: str) -> dict:
    """Parse ColabFold score files."""
    scores: dict[str, float] = {}

    # Try JSON scores file
    import json
    score_files = list(output_dir.glob(f"{job_name}*scores_rank_001*.json"))
    if not score_files:
        score_files = list(output_dir.glob(f"{job_name}*scores*.json"))

    if score_files:
        try:
            data = json.loads(score_files[0].read_text())
            scores["iptm"] = float(data.get("iptm", 0.0))
            scores["ptm"] = float(data.get("ptm", 0.0))
            plddt_list = data.get("plddt", [])
            if plddt_list:
                scores["plddt"] = sum(plddt_list) / len(plddt_list)
        except Exception as exc:
            logger.debug("Failed to parse ColabFold scores: %s", exc)

    return scores
