"""Tests for result browsing endpoints."""

from __future__ import annotations

import csv
import os
from pathlib import Path
from unittest.mock import patch

from service.models import Job


def _create_mock_results(job_id: str, output_base: str = "test_job_results") -> Path:
    """Create mock pipeline output files for testing."""
    job_dir = Path(output_base) / job_id / "20260301_run"
    scoring_dir = job_dir / "scoring"
    report_dir = job_dir / "report"
    scoring_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)

    # scored_candidates.csv
    csv_path = scoring_dir / "scored_candidates.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "candidate_id", "candidate_type", "modality", "source",
            "smiles", "sequence", "composite_score", "rank",
            "binding_score", "selectivity_score", "drug_likeness",
            "admet_score", "admet_flags", "delivery_system",
            "perturbation_score", "molecular_weight", "iupac_name",
        ])
        writer.writeheader()
        for i in range(5):
            writer.writerow({
                "candidate_id": f"CAND_{i:03d}",
                "candidate_type": "de_novo",
                "modality": "small_molecule",
                "source": "test",
                "smiles": f"C{'C' * i}O",
                "sequence": "",
                "composite_score": round(0.9 - i * 0.1, 2),
                "rank": i + 1,
                "binding_score": round(0.8 - i * 0.05, 2),
                "selectivity_score": 0.7,
                "drug_likeness": 0.6,
                "admet_score": 0.5,
                "admet_flags": 0,
                "delivery_system": "oral",
                "perturbation_score": 0.4,
                "molecular_weight": 300.0 + i * 10,
                "iupac_name": "",
            })

    # final_report.html
    (report_dir / "final_report.html").write_text("<html><body>Report</body></html>")

    return job_dir


def _submit_completed_job(client, auth_headers, test_db):
    """Submit a job and mark it completed with mock results."""
    with patch("service.routers.jobs_router.job_manager") as mock_jm:
        mock_jm.submit_job.return_value = True
        resp = client.post("/api/jobs/", json={
            "modality": "small_molecule",
            "mode": "antagonist",
            "target": "EGFR",
        }, headers=auth_headers)
    job_id = resp.json()["id"]

    result_dir = _create_mock_results(job_id)

    job = test_db.query(Job).filter(Job.id == job_id).first()
    job.status = "completed"
    job.result_dir = str(result_dir)
    test_db.commit()

    return job_id, result_dir


def test_get_result_summary(client, auth_headers, test_db):
    job_id, _ = _submit_completed_job(client, auth_headers, test_db)

    resp = client.get(f"/api/results/{job_id}", headers=auth_headers)
    assert resp.status_code == 200
    data = resp.json()
    assert data["job_id"] == job_id
    assert data["total_candidates"] == 5
    assert data["top_score"] == 0.9
    assert len(data["files"]) > 0


def test_get_candidates(client, auth_headers, test_db):
    job_id, _ = _submit_completed_job(client, auth_headers, test_db)

    resp = client.get(f"/api/results/{job_id}/candidates", headers=auth_headers)
    assert resp.status_code == 200
    data = resp.json()
    assert data["total"] == 5
    assert len(data["candidates"]) == 5
    assert data["candidates"][0]["composite_score"] == 0.9


def test_get_candidates_paginated(client, auth_headers, test_db):
    job_id, _ = _submit_completed_job(client, auth_headers, test_db)

    resp = client.get(
        f"/api/results/{job_id}/candidates?page=1&page_size=2",
        headers=auth_headers,
    )
    assert resp.status_code == 200
    data = resp.json()
    assert data["total"] == 5
    assert len(data["candidates"]) == 2


def test_result_not_found(client, auth_headers):
    resp = client.get("/api/results/nonexistent", headers=auth_headers)
    assert resp.status_code == 404


def test_download_file(client, auth_headers, test_db):
    job_id, result_dir = _submit_completed_job(client, auth_headers, test_db)

    resp = client.get(
        f"/api/results/{job_id}/files/report/final_report.html",
        headers=auth_headers,
    )
    assert resp.status_code == 200
    assert b"Report" in resp.content
