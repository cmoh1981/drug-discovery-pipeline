"""Tests for job submission and management endpoints."""

from __future__ import annotations

from unittest.mock import patch


def test_submit_job(client, auth_headers):
    with patch("service.routers.jobs_router.job_manager") as mock_jm:
        mock_jm.submit_job.return_value = True
        resp = client.post("/api/jobs/", json={
            "modality": "small_molecule",
            "mode": "antagonist",
            "target": "EGFR",
            "tissue": "lung",
            "top_n": 10,
        }, headers=auth_headers)
    assert resp.status_code == 202
    data = resp.json()
    assert data["target"] == "EGFR"
    assert data["modality"] == "small_molecule"
    assert data["status"] in ("pending", "running")
    assert "id" in data


def test_submit_job_invalid_modality(client, auth_headers):
    resp = client.post("/api/jobs/", json={
        "modality": "invalid",
        "mode": "antagonist",
        "target": "EGFR",
    }, headers=auth_headers)
    assert resp.status_code == 422


def test_submit_job_unauthorized(client):
    resp = client.post("/api/jobs/", json={
        "modality": "peptide",
        "mode": "agonist",
        "target": "TP53",
    })
    assert resp.status_code in (401, 403)


def test_list_jobs_empty(client, auth_headers):
    resp = client.get("/api/jobs/", headers=auth_headers)
    assert resp.status_code == 200
    data = resp.json()
    assert data["total"] == 0
    assert data["jobs"] == []


def test_list_jobs_with_data(client, auth_headers):
    with patch("service.routers.jobs_router.job_manager") as mock_jm:
        mock_jm.submit_job.return_value = True
        client.post("/api/jobs/", json={
            "modality": "peptide",
            "mode": "agonist",
            "target": "TP53",
        }, headers=auth_headers)
        client.post("/api/jobs/", json={
            "modality": "small_molecule",
            "mode": "antagonist",
            "target": "EGFR",
        }, headers=auth_headers)

    resp = client.get("/api/jobs/", headers=auth_headers)
    assert resp.status_code == 200
    data = resp.json()
    assert data["total"] == 2
    assert len(data["jobs"]) == 2


def test_get_job(client, auth_headers):
    with patch("service.routers.jobs_router.job_manager") as mock_jm:
        mock_jm.submit_job.return_value = True
        create_resp = client.post("/api/jobs/", json={
            "modality": "peptide",
            "mode": "agonist",
            "target": "BRCA1",
        }, headers=auth_headers)
    job_id = create_resp.json()["id"]

    resp = client.get(f"/api/jobs/{job_id}", headers=auth_headers)
    assert resp.status_code == 200
    assert resp.json()["target"] == "BRCA1"


def test_get_job_not_found(client, auth_headers):
    resp = client.get("/api/jobs/nonexistent-id", headers=auth_headers)
    assert resp.status_code == 404


def test_cancel_job(client, auth_headers):
    with patch("service.routers.jobs_router.job_manager") as mock_jm:
        mock_jm.submit_job.return_value = True
        mock_jm.cancel_job.return_value = True
        create_resp = client.post("/api/jobs/", json={
            "modality": "peptide",
            "mode": "agonist",
            "target": "MYC",
        }, headers=auth_headers)
    job_id = create_resp.json()["id"]

    with patch("service.routers.jobs_router.job_manager") as mock_jm:
        mock_jm.cancel_job.return_value = True
        resp = client.delete(f"/api/jobs/{job_id}", headers=auth_headers)
    assert resp.status_code == 204


def test_get_progress(client, auth_headers):
    with patch("service.routers.jobs_router.job_manager") as mock_jm:
        mock_jm.submit_job.return_value = True
        create_resp = client.post("/api/jobs/", json={
            "modality": "small_molecule",
            "mode": "antagonist",
            "target": "VEGFR",
        }, headers=auth_headers)
    job_id = create_resp.json()["id"]

    resp = client.get(f"/api/jobs/{job_id}/progress", headers=auth_headers)
    assert resp.status_code == 200
    data = resp.json()
    assert "step" in data
    assert "total" in data
    assert data["total"] == 10
