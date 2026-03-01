"""Cross-run comparison and config template endpoints."""

from __future__ import annotations

import json

from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.orm import Session

from service.comparison import compare_runs
from service.database import get_db
from service.dependencies import get_current_user
from service.models import ConfigTemplate, Job, User
from service.schemas import (
    CompareRequest,
    CompareResponse,
    ComparedCandidate,
    JobResponse,
    TemplateCreate,
    TemplateResponse,
)
from service.worker import job_manager

router = APIRouter(prefix="/api", tags=["compare"])


@router.post("/compare", response_model=CompareResponse)
def compare_jobs(
    body: CompareRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    # Verify all jobs belong to user
    result_dirs: dict[str, str | None] = {}
    for job_id in body.job_ids:
        job = db.query(Job).filter(Job.id == job_id, Job.user_id == current_user.id).first()
        if not job:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Job {job_id} not found",
            )
        if job.status != "completed":
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=f"Job {job_id} is not completed (status: {job.status})",
            )
        result_dirs[job_id] = job.result_dir

    result = compare_runs(body.job_ids, result_dirs, body.top_n)
    return CompareResponse(
        job_ids=result["job_ids"],
        shared_candidates=[ComparedCandidate(**c) for c in result["shared_candidates"]],
        unique_per_job=result["unique_per_job"],
        total_compared=result["total_compared"],
    )


@router.post("/jobs/{job_id}/rerun", response_model=JobResponse, status_code=status.HTTP_202_ACCEPTED)
def rerun_job(
    job_id: str,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    original = db.query(Job).filter(Job.id == job_id, Job.user_id == current_user.id).first()
    if not original:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Job not found")

    new_job = Job(
        user_id=current_user.id,
        modality=original.modality,
        mode=original.mode,
        target=original.target,
        tissue=original.tissue,
        top_n=original.top_n,
        use_runpod=original.use_runpod,
        config_overrides=original.config_overrides,
    )
    db.add(new_job)
    db.commit()
    db.refresh(new_job)

    job_manager.submit_job(new_job.id, {
        "modality": new_job.modality,
        "mode": new_job.mode,
        "target": new_job.target,
    })
    return new_job


@router.get("/jobs/templates", response_model=list[TemplateResponse])
def list_templates(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    return db.query(ConfigTemplate).filter(ConfigTemplate.user_id == current_user.id).all()


@router.post("/jobs/templates", response_model=TemplateResponse, status_code=status.HTTP_201_CREATED)
def create_template(
    body: TemplateCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    template = ConfigTemplate(
        user_id=current_user.id,
        name=body.name,
        description=body.description,
        config_json=json.dumps(body.config),
    )
    db.add(template)
    db.commit()
    db.refresh(template)
    return template
