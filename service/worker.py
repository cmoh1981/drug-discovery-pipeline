"""Background job execution wrapping the existing pipeline."""

from __future__ import annotations

import json
import logging
import traceback
from concurrent.futures import Future, ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from threading import Lock

from service.config import settings

logger = logging.getLogger(__name__)


class JobManager:
    """Manages pipeline job execution in a thread pool."""

    def __init__(self, max_workers: int | None = None):
        self._max_workers = max_workers or settings.MAX_CONCURRENT_JOBS
        self._executor = ThreadPoolExecutor(max_workers=self._max_workers)
        self._futures: dict[str, Future] = {}
        self._lock = Lock()

    @property
    def active_count(self) -> int:
        with self._lock:
            return sum(1 for f in self._futures.values() if not f.done())

    def submit_job(self, job_id: str, job_params: dict) -> bool:
        """Submit a pipeline job for background execution.

        Returns True if submitted, False if at capacity.
        """
        if self.active_count >= self._max_workers:
            return False

        future = self._executor.submit(self._run_job, job_id, job_params)
        with self._lock:
            self._futures[job_id] = future
        return True

    def cancel_job(self, job_id: str) -> bool:
        with self._lock:
            future = self._futures.get(job_id)
        if future and not future.done():
            return future.cancel()
        return False

    def shutdown(self, wait: bool = True) -> None:
        self._executor.shutdown(wait=wait)

    def _run_job(self, job_id: str, job_params: dict) -> Path:
        """Execute the pipeline in a worker thread.

        Updates the job record in DB via direct session (thread-safe).
        """
        from sqlalchemy.orm import Session
        from service.database import SessionLocal
        from service.models import Job

        db: Session = SessionLocal()
        try:
            job = db.query(Job).filter(Job.id == job_id).first()
            if not job:
                raise ValueError(f"Job {job_id} not found")

            # Mark running
            job.status = "running"
            job.started_at = datetime.now(timezone.utc)
            db.commit()

            # Build PipelineConfig from job params
            from drugdiscovery.types import Modality, ModeOfAction, PipelineConfig

            output_dir = str(settings.output_path / job_id)
            Path(output_dir).mkdir(parents=True, exist_ok=True)

            cfg = PipelineConfig(
                modality=Modality(job.modality),
                mode=ModeOfAction(job.mode),
                target=job.target,
                tissue=job.tissue or "",
                top_n=job.top_n,
                use_runpod=job.use_runpod,
                output_dir=output_dir,
            )

            # Apply config overrides
            if job.config_overrides:
                overrides = json.loads(job.config_overrides)
                if "scoring_weights" in overrides and isinstance(overrides["scoring_weights"], dict):
                    cfg.scoring_weights.update(overrides["scoring_weights"])
                if "top_n" in overrides:
                    cfg.top_n = int(overrides["top_n"])

            # Run the pipeline (existing function, zero changes)
            from drugdiscovery.pipeline import run_pipeline

            result_path = run_pipeline(cfg)

            # Mark completed
            job = db.query(Job).filter(Job.id == job_id).first()
            job.status = "completed"
            job.completed_at = datetime.now(timezone.utc)
            job.result_dir = str(result_path)

            # Extract summary stats from scored_candidates CSV if available
            self._update_job_stats(job, result_path)

            db.commit()
            logger.info("Job %s completed: %s", job_id, result_path)
            return result_path

        except Exception as exc:
            logger.error("Job %s failed: %s", job_id, exc, exc_info=True)
            db.rollback()
            job = db.query(Job).filter(Job.id == job_id).first()
            if job:
                job.status = "failed"
                job.error_message = f"{type(exc).__name__}: {exc}"
                job.completed_at = datetime.now(timezone.utc)
                db.commit()
            raise
        finally:
            db.close()

    @staticmethod
    def _update_job_stats(job, result_path: Path) -> None:
        """Read scored_candidates.csv to populate summary stats."""
        import csv

        csv_candidates = list(result_path.rglob("scored_candidates.csv"))
        if not csv_candidates:
            return
        csv_path = csv_candidates[0]
        try:
            with open(csv_path, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                rows = list(reader)
            if rows:
                job.total_candidates = len(rows)
                scores = [float(r.get("composite_score", 0)) for r in rows]
                job.top_score = max(scores) if scores else None
        except Exception:
            pass


# Singleton instance
job_manager = JobManager()
