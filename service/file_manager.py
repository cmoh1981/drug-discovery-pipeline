"""Read pipeline output files as structured data."""

from __future__ import annotations

import csv
from pathlib import Path

from service.config import settings


def get_job_output_dir(job_id: str, result_dir: str | None = None) -> Path | None:
    """Resolve the output directory for a completed job."""
    if result_dir:
        p = Path(result_dir)
        if p.is_dir():
            return p
    # Fallback: scan job_results/{job_id}
    job_dir = settings.output_path / job_id
    if job_dir.is_dir():
        # The pipeline creates a timestamped subdirectory
        subdirs = sorted(job_dir.iterdir())
        if subdirs:
            return subdirs[0] if subdirs[0].is_dir() else job_dir
    return None


def list_job_files(output_dir: Path) -> list[str]:
    """List all files in a job's output directory (relative paths)."""
    if not output_dir.is_dir():
        return []
    return sorted(
        str(f.relative_to(output_dir))
        for f in output_dir.rglob("*")
        if f.is_file()
    )


def get_candidates_data(
    output_dir: Path,
    page: int = 1,
    page_size: int = 20,
    sort_by: str = "composite_score",
    descending: bool = True,
) -> tuple[list[dict], int]:
    """Read scored_candidates.csv and return paginated results.

    Returns (rows, total_count).
    """
    csv_files = list(output_dir.rglob("scored_candidates.csv"))
    if not csv_files:
        return [], 0

    csv_path = csv_files[0]
    rows = _read_csv(csv_path)

    # Sort
    if rows and sort_by in rows[0]:
        try:
            rows.sort(key=lambda r: float(r.get(sort_by, 0) or 0), reverse=descending)
        except (ValueError, TypeError):
            rows.sort(key=lambda r: r.get(sort_by, ""), reverse=descending)

    total = len(rows)
    start = (page - 1) * page_size
    end = start + page_size
    return rows[start:end], total


def get_job_summary(output_dir: Path) -> dict:
    """Build a summary dict from pipeline output files."""
    files = list_job_files(output_dir)
    csv_files = list(output_dir.rglob("scored_candidates.csv"))
    total_candidates = 0
    top_score = None

    if csv_files:
        rows = _read_csv(csv_files[0])
        total_candidates = len(rows)
        scores = [float(r.get("composite_score", 0) or 0) for r in rows]
        if scores:
            top_score = max(scores)

    return {
        "total_candidates": total_candidates,
        "top_score": top_score,
        "files": files,
    }


def _read_csv(path: Path) -> list[dict]:
    """Read a CSV into a list of dicts."""
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f))
