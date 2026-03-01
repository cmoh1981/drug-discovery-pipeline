"""Cross-run candidate comparison logic."""

from __future__ import annotations

from collections import defaultdict

from service.file_manager import get_candidates_data, get_job_output_dir


def compare_runs(
    job_ids: list[str],
    result_dirs: dict[str, str | None],
    top_n: int = 20,
) -> dict:
    """Compare candidates across multiple pipeline runs.

    Matches candidates by SMILES (small molecules) or sequence (peptides).
    Returns shared candidates with per-run scores and unique counts.
    """
    # Load candidates per job
    job_candidates: dict[str, list[dict]] = {}
    for job_id in job_ids:
        output_dir = get_job_output_dir(job_id, result_dirs.get(job_id))
        if not output_dir:
            job_candidates[job_id] = []
            continue
        rows, _ = get_candidates_data(output_dir, page=1, page_size=top_n)
        job_candidates[job_id] = rows

    # Build identifier → {job_id: score} mapping
    candidate_map: dict[str, dict[str, float | None]] = defaultdict(dict)
    job_identifiers: dict[str, set[str]] = {jid: set() for jid in job_ids}

    for job_id, candidates in job_candidates.items():
        for c in candidates:
            identifier = _get_candidate_identifier(c)
            if not identifier:
                continue
            score = _safe_float(c.get("composite_score"))
            candidate_map[identifier][job_id] = score
            job_identifiers[job_id].add(identifier)

    # Find shared candidates (appear in 2+ runs)
    shared = []
    for identifier, scores in candidate_map.items():
        if len(scores) >= 2:
            # Pad missing jobs with None
            full_scores = {jid: scores.get(jid) for jid in job_ids}
            shared.append({
                "identifier": identifier,
                "scores": full_scores,
            })

    # Sort shared by average score (descending)
    shared.sort(
        key=lambda x: _avg_score(x["scores"]),
        reverse=True,
    )

    # Count unique per job
    all_shared_ids = {s["identifier"] for s in shared}
    unique_per_job = {
        jid: len(ids - all_shared_ids)
        for jid, ids in job_identifiers.items()
    }

    return {
        "job_ids": job_ids,
        "shared_candidates": shared[:top_n],
        "unique_per_job": unique_per_job,
        "total_compared": sum(len(v) for v in job_candidates.values()),
    }


def _get_candidate_identifier(candidate: dict) -> str:
    """Get the best identifier for matching: SMILES for SM, sequence for peptides."""
    smiles = candidate.get("smiles", "").strip()
    if smiles:
        return f"SMILES:{smiles}"
    sequence = candidate.get("sequence", "").strip()
    if sequence:
        return f"SEQ:{sequence}"
    return candidate.get("candidate_id", "")


def _safe_float(value) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (ValueError, TypeError):
        return None


def _avg_score(scores: dict[str, float | None]) -> float:
    valid = [v for v in scores.values() if v is not None]
    return sum(valid) / len(valid) if valid else 0.0
