"""Result browsing and log streaming endpoints."""

from __future__ import annotations

from fastapi import APIRouter, Depends, HTTPException, Query, status
from fastapi.responses import FileResponse, RedirectResponse, StreamingResponse
from sqlalchemy.orm import Session

from service.database import get_db
from service.dependencies import get_current_user
from service.file_manager import get_candidates_data, get_job_output_dir, get_job_summary
from service.log_streamer import stream_logs
from service.models import Job, User
from service.schemas import CandidateListResponse, CandidateResponse, JobSummary

router = APIRouter(prefix="/api/results", tags=["results"])


def _get_user_job(job_id: str, db: Session, current_user: User) -> Job:
    job = db.query(Job).filter(Job.id == job_id, Job.user_id == current_user.id).first()
    if not job:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Job not found")
    return job


@router.get("/{job_id}", response_model=JobSummary)
def get_result_summary(
    job_id: str,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    job = _get_user_job(job_id, db, current_user)
    output_dir = get_job_output_dir(job_id, job.result_dir)
    if not output_dir:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Result directory not found",
        )
    summary = get_job_summary(output_dir)
    return JobSummary(
        job_id=job.id,
        status=job.status,
        target=job.target,
        modality=job.modality,
        mode=job.mode,
        total_candidates=summary["total_candidates"],
        top_score=summary["top_score"],
        files=summary["files"],
        started_at=job.started_at,
        completed_at=job.completed_at,
    )


@router.get("/{job_id}/candidates", response_model=CandidateListResponse)
def get_candidates(
    job_id: str,
    page: int = Query(1, ge=1),
    page_size: int = Query(20, ge=1, le=100),
    sort_by: str = Query("composite_score"),
    descending: bool = Query(True),
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    job = _get_user_job(job_id, db, current_user)
    output_dir = get_job_output_dir(job_id, job.result_dir)
    if not output_dir:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Result directory not found",
        )
    rows, total = get_candidates_data(output_dir, page, page_size, sort_by, descending)
    candidates = [CandidateResponse(**r) for r in rows]
    return CandidateListResponse(
        candidates=candidates, total=total, page=page, page_size=page_size,
    )


@router.get("/{job_id}/dashboard")
def get_dashboard(
    job_id: str,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    job = _get_user_job(job_id, db, current_user)
    output_dir = get_job_output_dir(job_id, job.result_dir)
    if not output_dir:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Result directory not found",
        )
    html_files = list(output_dir.rglob("final_report.html"))
    if not html_files:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Dashboard HTML not found",
        )
    # Redirect to the static mount
    relative = html_files[0].relative_to(output_dir)
    return RedirectResponse(url=f"/reports/{job_id}/{relative}")


@router.get("/{job_id}/files/{file_path:path}")
def download_file(
    job_id: str,
    file_path: str,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    job = _get_user_job(job_id, db, current_user)
    output_dir = get_job_output_dir(job_id, job.result_dir)
    if not output_dir:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Result directory not found",
        )
    target_file = output_dir / file_path
    # Prevent path traversal
    try:
        target_file.resolve().relative_to(output_dir.resolve())
    except ValueError:
        raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Access denied")
    if not target_file.is_file():
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="File not found")
    return FileResponse(target_file)


@router.get("/{job_id}/logs")
def stream_job_logs(
    job_id: str,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    _get_user_job(job_id, db, current_user)
    return StreamingResponse(
        stream_logs(job_id),
        media_type="text/event-stream",
        headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"},
    )
