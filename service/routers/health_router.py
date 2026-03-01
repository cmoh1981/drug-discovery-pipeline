"""Health check endpoint."""

from __future__ import annotations

from fastapi import APIRouter, Depends
from sqlalchemy import text
from sqlalchemy.orm import Session

from service.database import get_db
from service.schemas import HealthResponse
from service.worker import job_manager

router = APIRouter(tags=["health"])


@router.get("/api/health", response_model=HealthResponse)
def health_check(db: Session = Depends(get_db)):
    db_status = "connected"
    try:
        db.execute(text("SELECT 1"))
    except Exception:
        db_status = "disconnected"
    return HealthResponse(
        status="ok" if db_status == "connected" else "degraded",
        database=db_status,
        active_jobs=job_manager.active_count,
    )
