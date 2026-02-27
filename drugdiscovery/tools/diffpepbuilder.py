"""DiffPepBuilder wrapper for structure-based peptide design.

Refactored from YARS2 pipeline: scripts/run_diffpepbuilder.py

Supports two execution backends:
  - Local subprocess: DiffPepBuilder cloned at tools/DiffPepBuilder/ (or a
    caller-supplied path stored in the environment variable
    DIFFPEPBUILDER_DIR).
  - RunPod serverless: GPU endpoint via drugdiscovery.utils.compute.RunPodClient.

PDB coordinate utilities delegate to drugdiscovery.utils.io so that the
three-to-one residue map and CA-extraction logic stay in a single place.
"""

from __future__ import annotations

import json
import logging
import os
import sys
from pathlib import Path
from typing import Optional

from drugdiscovery.utils.io import extract_sequence_from_pdb, parse_pdb_binding_coords

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_DEFAULT_DPB_DIR = Path(os.environ.get("DIFFPEPBUILDER_DIR", "tools/DiffPepBuilder"))

# Entry-point scripts tried in order when locating DiffPepBuilder.
_DPB_ENTRY_SCRIPTS = [
    "inference.py",
    "run_inference.py",
    "sample.py",
    "generate.py",
    "main.py",
]

# B-factor column slice in a PDB ATOM record (columns 61-66, 0-indexed 60:66).
_BFACTOR_SLICE = slice(60, 66)


# ---------------------------------------------------------------------------
# Internal PDB helpers
# ---------------------------------------------------------------------------

def _read_pdb_text(pdb_path: str) -> str:
    """Read a PDB file and return its full text content."""
    return Path(pdb_path).read_text(encoding="utf-8")


def _mean_bfactor(pdb_path: str) -> Optional[float]:
    """Compute mean B-factor from ATOM records (proxy for pLDDT in AF models)."""
    bfactors: list[float] = []
    for line in _read_pdb_text(pdb_path).splitlines():
        if not line.startswith("ATOM"):
            continue
        try:
            bfactors.append(float(line[_BFACTOR_SLICE].strip()))
        except (ValueError, IndexError):
            continue
    return sum(bfactors) / len(bfactors) if bfactors else None


def extract_sequence_from_pdb_extended(pdb_path: str, chain: str = "B") -> str:
    """Extract the peptide amino acid sequence from a generated PDB file.

    Tries SEQRES records first; falls back to ATOM records for *chain*.
    This function overrides the import from utils.io with an augmented version
    that also checks SEQRES records, which DiffPepBuilder output often contains.

    Args:
        pdb_path: Absolute path to the PDB file.
        chain:    Chain identifier for the peptide (default "B").

    Returns:
        One-letter amino acid sequence string.  Unknown residues become "X".
    """
    three_to_one = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    }
    pdb_text = _read_pdb_text(pdb_path)

    # --- SEQRES records (preferred: explicit sequence annotation) ---
    seqres_parts: list[str] = []
    for line in pdb_text.splitlines():
        if line.startswith("SEQRES"):
            record_chain = line[11] if len(line) > 11 else ""
            if record_chain != chain and chain != "":
                continue
            for token in line[19:].split():
                seqres_parts.append(three_to_one.get(token.upper(), "X"))
    if seqres_parts:
        return "".join(seqres_parts)

    # --- ATOM records fallback (delegate to utils.io implementation) ---
    from drugdiscovery.utils.io import extract_sequence_from_pdb as _util_extract  # noqa: PLC0415
    return _util_extract(pdb_text, chain=chain)


# ---------------------------------------------------------------------------
# Result extraction from raw DiffPepBuilder output
# ---------------------------------------------------------------------------

def _extract_results_from_dir(raw_output_dir: Path) -> list[dict]:
    """Parse generated files from a DiffPepBuilder local output directory.

    Tries, in order:
      1. results.csv  — preferred structured output
      2. sequences.txt — plain sequence list
      3. *.pdb files  — extract sequence + B-factor confidence

    Returns list of dicts: {sequence, pdb_path, confidence}.
    """
    candidates: list[dict] = []
    counter = 0

    # 1. Structured CSV
    results_csv = raw_output_dir / "results.csv"
    if results_csv.exists():
        import csv  # noqa: PLC0415
        logger.info("[DiffPepBuilder] Reading results from %s", results_csv)
        with open(results_csv, newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                counter += 1
                seq = row.get("sequence", row.get("peptide_sequence", "")).strip().upper()
                try:
                    conf = float(row.get("confidence", row.get("score", "nan")))
                except (ValueError, TypeError):
                    conf = float("nan")
                candidates.append({
                    "sequence": seq,
                    "pdb_path": row.get("pdb_path", ""),
                    "confidence": conf,
                })
        return candidates

    # 2. Plain sequence list
    seq_file = raw_output_dir / "sequences.txt"
    if seq_file.exists():
        logger.info("[DiffPepBuilder] Reading sequences from %s", seq_file)
        for line in seq_file.read_text(encoding="utf-8").splitlines():
            seq = line.strip().upper()
            if not seq or seq.startswith("#"):
                continue
            counter += 1
            candidates.append({"sequence": seq, "pdb_path": "", "confidence": float("nan")})
        return candidates

    # 3. Individual PDB files
    pdb_files = sorted(raw_output_dir.glob("*.pdb"))
    logger.info(
        "[DiffPepBuilder] Parsing %d PDB file(s) from %s", len(pdb_files), raw_output_dir
    )
    for pdb_file in pdb_files:
        seq = extract_sequence_from_pdb_extended(str(pdb_file), chain="B")
        conf_raw = _mean_bfactor(str(pdb_file))
        conf = conf_raw if conf_raw is not None else float("nan")
        candidates.append({
            "sequence": seq,
            "pdb_path": str(pdb_file),
            "confidence": conf,
        })

    logger.info(
        "[DiffPepBuilder] Extracted %d candidate(s) from output directory.",
        len(candidates),
    )
    return candidates


# ---------------------------------------------------------------------------
# Local subprocess execution
# ---------------------------------------------------------------------------

def _find_entry_script(dpb_dir: Path) -> Path:
    """Locate the DiffPepBuilder inference entry-point script."""
    for name in _DPB_ENTRY_SCRIPTS:
        candidate = dpb_dir / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        f"[DiffPepBuilder] Could not find an entry-point script in '{dpb_dir}'. "
        f"Tried: {_DPB_ENTRY_SCRIPTS}. "
        "Clone DiffPepBuilder with: "
        "git clone https://github.com/YuzheWangPKU/DiffPepBuilder tools/DiffPepBuilder"
    )


def run_diffpepbuilder_local(
    target_pdb: str,
    binding_site_residues: list[int],
    output_dir: str,
    num_samples: int = 100,
    peptide_length_min: int = 8,
    peptide_length_max: int = 25,
    dpb_dir: Optional[str] = None,
) -> list[dict]:
    """Run DiffPepBuilder locally via subprocess.

    Writes a JSON config file to *output_dir*, then invokes DiffPepBuilder's
    inference script as a child process.  The raw output directory is parsed
    and a clean candidate list is returned.

    Args:
        target_pdb:            Path to the target protein PDB file.
        binding_site_residues: List of residue numbers defining the binding site.
        output_dir:            Directory for all output files.
        num_samples:           Number of peptide structures to sample.
        peptide_length_min:    Minimum peptide length in residues.
        peptide_length_max:    Maximum peptide length in residues.
        dpb_dir:               Path to the cloned DiffPepBuilder repository.
                               Falls back to the DIFFPEPBUILDER_DIR environment
                               variable and then to "tools/DiffPepBuilder/".

    Returns:
        List of dicts with keys:
          ``sequence``   – amino acid sequence (str)
          ``pdb_path``   – path to generated PDB file, empty string if unavailable
          ``confidence`` – mean B-factor / pLDDT proxy (float, may be NaN)
    """
    import subprocess  # noqa: PLC0415

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    diffpepbuilder_dir = Path(dpb_dir) if dpb_dir else _DEFAULT_DPB_DIR
    entry_script = _find_entry_script(diffpepbuilder_dir)

    # Extract binding site CA coordinates for the config
    pdb_text = Path(target_pdb).read_text(encoding="utf-8")
    binding_coords = parse_pdb_binding_coords(pdb_text, binding_site_residues, chain="A")
    if not binding_coords:
        raise ValueError(
            f"[DiffPepBuilder] No CA atoms found for residues {binding_site_residues} "
            f"in '{target_pdb}'. Check chain ID and residue numbering."
        )

    # Compute centroid
    n = len(binding_coords)
    centroid = {
        "x": round(sum(c["x"] for c in binding_coords) / n, 3),
        "y": round(sum(c["y"] for c in binding_coords) / n, 3),
        "z": round(sum(c["z"] for c in binding_coords) / n, 3),
    }

    config = {
        "target_pdb": str(Path(target_pdb).resolve()),
        "binding_site_residues": binding_site_residues,
        "binding_site_coords": binding_coords,
        "binding_site_centroid": centroid,
        "peptide_length_min": peptide_length_min,
        "peptide_length_max": peptide_length_max,
        "num_samples": num_samples,
        "output_dir": str(out_dir.resolve()),
    }

    config_path = out_dir / "diffpepbuilder_input_config.json"
    config_path.write_text(json.dumps(config, indent=2), encoding="utf-8")
    logger.info("[DiffPepBuilder] Config written to %s", config_path)

    raw_output_dir = out_dir / "diffpepbuilder_raw_output"
    raw_output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        sys.executable,
        str(entry_script),
        "--config", str(config_path),
        "--target_pdb", config["target_pdb"],
        "--hotspots", ",".join(str(r) for r in binding_site_residues),
        "--pep_length_min", str(peptide_length_min),
        "--pep_length_max", str(peptide_length_max),
        "--num_samples", str(num_samples),
        "--output_dir", str(raw_output_dir),
    ]

    logger.info("[DiffPepBuilder] Launching subprocess: %s", " ".join(cmd))
    result = subprocess.run(
        cmd,
        cwd=str(diffpepbuilder_dir),
        capture_output=False,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"[DiffPepBuilder] Subprocess exited with code {result.returncode}."
        )

    logger.info("[DiffPepBuilder] Subprocess completed; parsing output.")
    candidates = _extract_results_from_dir(raw_output_dir)
    logger.info("[DiffPepBuilder] Local run produced %d candidate(s).", len(candidates))
    return candidates


# ---------------------------------------------------------------------------
# RunPod serverless execution
# ---------------------------------------------------------------------------

def run_diffpepbuilder_runpod(
    target_pdb: str,
    binding_site_residues: list[int],
    output_dir: str,
    num_samples: int = 100,
    peptide_length_min: int = 8,
    peptide_length_max: int = 25,
    endpoint_id: str = "",
) -> list[dict]:
    """Run DiffPepBuilder via a RunPod serverless GPU endpoint.

    Reads PDB coordinates and centroid locally, then submits the configuration
    payload to the RunPod endpoint.  The endpoint is expected to return a JSON
    object with the keys ``sequences``, ``pdbs`` (optional), and
    ``confidences`` (optional).

    Args:
        target_pdb:            Path to the target protein PDB file.
        binding_site_residues: List of residue numbers defining the binding site.
        output_dir:            Directory to write downloaded PDB files.
        num_samples:           Number of peptide structures to sample.
        peptide_length_min:    Minimum peptide length in residues.
        peptide_length_max:    Maximum peptide length in residues.
        endpoint_id:           RunPod serverless endpoint ID.  Falls back to the
                               RUNPOD_ENDPOINT_ID environment variable.

    Returns:
        List of dicts with keys:
          ``sequence``   – amino acid sequence (str)
          ``pdb_path``   – path to the saved PDB file, empty if not returned
          ``confidence`` – model confidence score (float, may be NaN)

    Raises:
        RuntimeError: If the RunPod job fails or RUNPOD_API_KEY is not set.
        ValueError:   If no binding site coordinates can be extracted.
    """
    from drugdiscovery.utils.compute import RunPodClient  # noqa: PLC0415

    resolved_endpoint = endpoint_id or os.environ.get("RUNPOD_ENDPOINT_ID", "")
    if not resolved_endpoint:
        raise RuntimeError(
            "[DiffPepBuilder] RunPod endpoint ID is required. "
            "Pass endpoint_id or set the RUNPOD_ENDPOINT_ID environment variable."
        )

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pdb_text = Path(target_pdb).read_text(encoding="utf-8")
    binding_coords = parse_pdb_binding_coords(pdb_text, binding_site_residues, chain="A")
    if not binding_coords:
        raise ValueError(
            f"[DiffPepBuilder] No CA atoms found for residues {binding_site_residues} "
            f"in '{target_pdb}'."
        )

    n = len(binding_coords)
    centroid = {
        "x": round(sum(c["x"] for c in binding_coords) / n, 3),
        "y": round(sum(c["y"] for c in binding_coords) / n, 3),
        "z": round(sum(c["z"] for c in binding_coords) / n, 3),
    }

    payload = {
        "target_pdb_content": pdb_text,
        "binding_site_residues": binding_site_residues,
        "binding_site_coords": binding_coords,
        "binding_site_centroid": centroid,
        "peptide_length_min": peptide_length_min,
        "peptide_length_max": peptide_length_max,
        "num_samples": num_samples,
    }

    logger.info(
        "[DiffPepBuilder] Submitting RunPod job to endpoint '%s' "
        "(%d hotspots, %d samples).",
        resolved_endpoint,
        len(binding_site_residues),
        num_samples,
    )

    client = RunPodClient()
    output = client.run_job(endpoint_id=resolved_endpoint, payload=payload)

    sequences: list[str] = output.get("sequences", [])
    pdbs: list[str] = output.get("pdbs", [])
    confidences: list = output.get("confidences", [])

    candidates: list[dict] = []
    for i, seq in enumerate(sequences):
        pdb_content = pdbs[i] if i < len(pdbs) else ""
        raw_conf = confidences[i] if i < len(confidences) else None

        try:
            conf = float(raw_conf) if raw_conf is not None else float("nan")
        except (ValueError, TypeError):
            conf = float("nan")

        saved_pdb = ""
        if pdb_content:
            pdb_file = out_dir / f"runpod_sample_{i + 1:05d}.pdb"
            pdb_file.write_text(pdb_content, encoding="utf-8")
            saved_pdb = str(pdb_file)

        candidates.append({
            "sequence": seq.strip().upper(),
            "pdb_path": saved_pdb,
            "confidence": conf,
        })

    logger.info(
        "[DiffPepBuilder] RunPod job returned %d candidate(s).", len(candidates)
    )
    return candidates
