"""Tests for cross-run comparison logic."""

from __future__ import annotations

import csv
from pathlib import Path
from unittest.mock import patch

from service.comparison import compare_runs, _get_candidate_identifier
from service.models import Job


def _create_run_results(job_id: str, candidates: list[dict], base: str = "test_job_results") -> str:
    """Create mock scored_candidates.csv for a job."""
    job_dir = Path(base) / job_id / "20260301_run" / "scoring"
    job_dir.mkdir(parents=True, exist_ok=True)
    csv_path = job_dir / "scored_candidates.csv"

    if candidates:
        fieldnames = list(candidates[0].keys())
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(candidates)

    return str(job_dir.parent)


def test_get_candidate_identifier_smiles():
    assert _get_candidate_identifier({"smiles": "CCO", "sequence": ""}) == "SMILES:CCO"


def test_get_candidate_identifier_sequence():
    assert _get_candidate_identifier({"smiles": "", "sequence": "ACGT"}) == "SEQ:ACGT"


def test_compare_shared_candidates():
    c1 = [
        {"candidate_id": "A", "smiles": "CCO", "sequence": "", "composite_score": "0.9"},
        {"candidate_id": "B", "smiles": "CCCO", "sequence": "", "composite_score": "0.8"},
    ]
    c2 = [
        {"candidate_id": "C", "smiles": "CCO", "sequence": "", "composite_score": "0.85"},
        {"candidate_id": "D", "smiles": "CCCCO", "sequence": "", "composite_score": "0.7"},
    ]

    dir1 = _create_run_results("job1", c1)
    dir2 = _create_run_results("job2", c2)

    result = compare_runs(
        ["job1", "job2"],
        {"job1": dir1, "job2": dir2},
        top_n=10,
    )
    assert len(result["shared_candidates"]) == 1
    assert result["shared_candidates"][0]["identifier"] == "SMILES:CCO"
    assert result["unique_per_job"]["job1"] == 1
    assert result["unique_per_job"]["job2"] == 1


def test_compare_no_shared():
    c1 = [{"candidate_id": "A", "smiles": "CCO", "sequence": "", "composite_score": "0.9"}]
    c2 = [{"candidate_id": "B", "smiles": "CCCO", "sequence": "", "composite_score": "0.8"}]

    dir1 = _create_run_results("job3", c1)
    dir2 = _create_run_results("job4", c2)

    result = compare_runs(
        ["job3", "job4"],
        {"job3": dir1, "job4": dir2},
        top_n=10,
    )
    assert len(result["shared_candidates"]) == 0


def test_compare_endpoint(client, auth_headers, test_db):
    """Test the /api/compare endpoint."""
    # Create two completed jobs
    with patch("service.routers.jobs_router.job_manager") as mock_jm:
        mock_jm.submit_job.return_value = True
        resp1 = client.post("/api/jobs/", json={
            "modality": "small_molecule",
            "mode": "antagonist",
            "target": "EGFR",
        }, headers=auth_headers)
        resp2 = client.post("/api/jobs/", json={
            "modality": "small_molecule",
            "mode": "antagonist",
            "target": "EGFR",
        }, headers=auth_headers)

    job1_id = resp1.json()["id"]
    job2_id = resp2.json()["id"]

    # Create mock results
    candidates = [
        {"candidate_id": "X", "smiles": "CCO", "sequence": "", "composite_score": "0.9"},
    ]
    dir1 = _create_run_results(job1_id, candidates)
    dir2 = _create_run_results(job2_id, candidates)

    # Mark completed
    for jid, rdir in [(job1_id, dir1), (job2_id, dir2)]:
        job = test_db.query(Job).filter(Job.id == jid).first()
        job.status = "completed"
        job.result_dir = rdir
    test_db.commit()

    resp = client.post("/api/compare", json={
        "job_ids": [job1_id, job2_id],
        "top_n": 10,
    }, headers=auth_headers)
    assert resp.status_code == 200
    data = resp.json()
    assert len(data["shared_candidates"]) == 1
