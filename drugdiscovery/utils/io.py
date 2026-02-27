"""Unified I/O utilities for CSV, FASTA, JSON, and PDB files."""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from typing import Any, Sequence

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# CSV
# ---------------------------------------------------------------------------

def read_csv(path: str | Path) -> list[dict[str, str]]:
    """Read a CSV file and return list of row dicts."""
    with open(path, "r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: str | Path, rows: Sequence[dict], fieldnames: list[str] | None = None) -> None:
    """Write rows to CSV. Infers fieldnames from first row if not given."""
    if not rows:
        logger.warning("write_csv: no rows to write to %s", path)
        return
    if fieldnames is None:
        fieldnames = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    logger.info("Wrote %d rows to %s", len(rows), path)


# ---------------------------------------------------------------------------
# FASTA
# ---------------------------------------------------------------------------

def read_fasta(path: str | Path) -> list[tuple[str, str]]:
    """Read FASTA file. Returns list of (header, sequence) tuples."""
    entries: list[tuple[str, str]] = []
    header = ""
    seq_parts: list[str] = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    entries.append((header, "".join(seq_parts)))
                header = line[1:].strip()
                seq_parts = []
            elif line:
                seq_parts.append(line)
    if header:
        entries.append((header, "".join(seq_parts)))
    return entries


def write_fasta(
    path: str | Path,
    entries: Sequence[tuple[str, str]],
    line_width: int = 80,
) -> None:
    """Write sequences to FASTA format."""
    with open(path, "w", encoding="utf-8") as fh:
        for header, seq in entries:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), line_width):
                fh.write(seq[i : i + line_width] + "\n")
    logger.info("Wrote %d sequences to %s", len(entries), path)


# ---------------------------------------------------------------------------
# JSON
# ---------------------------------------------------------------------------

def read_json(path: str | Path) -> Any:
    """Read a JSON file."""
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def write_json(path: str | Path, data: Any, indent: int = 2) -> None:
    """Write data to JSON file."""
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=indent, default=str)
    logger.info("Wrote JSON to %s", path)


# ---------------------------------------------------------------------------
# PDB helpers
# ---------------------------------------------------------------------------

def extract_sequence_from_pdb(pdb_text: str, chain: str = "A") -> str:
    """Extract amino acid sequence from PDB ATOM records for a given chain."""
    three_to_one = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    }
    seen_residues: dict[int, str] = {}
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[21] == chain:
            resnum = int(line[22:26].strip())
            resname = line[17:20].strip()
            if resnum not in seen_residues:
                aa = three_to_one.get(resname, "X")
                seen_residues[resnum] = aa
    return "".join(seen_residues[k] for k in sorted(seen_residues))


def parse_pdb_binding_coords(
    pdb_text: str,
    residue_numbers: list[int],
    chain: str = "A",
) -> list[dict]:
    """Extract CA atom coordinates for specified residues."""
    coords = []
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[21] != chain:
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        resnum = int(line[22:26].strip())
        if resnum in residue_numbers:
            coords.append({
                "resnum": resnum,
                "resname": line[17:20].strip(),
                "chain": chain,
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
            })
    return coords


def compute_centroid(coords: list[dict]) -> tuple[float, float, float]:
    """Compute centroid from a list of coordinate dicts."""
    if not coords:
        return (0.0, 0.0, 0.0)
    n = len(coords)
    cx = sum(c["x"] for c in coords) / n
    cy = sum(c["y"] for c in coords) / n
    cz = sum(c["z"] for c in coords) / n
    return (round(cx, 3), round(cy, 3), round(cz, 3))
