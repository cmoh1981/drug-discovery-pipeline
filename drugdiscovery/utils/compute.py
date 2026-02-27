"""GPU detection and RunPod cloud compute utilities.

Refactored from YARS2 pipeline: runpod_utils.py
"""

from __future__ import annotations

import logging
import os
import time
from typing import Any, Optional

logger = logging.getLogger(__name__)


def detect_device(requested: str = "auto") -> str:
    """Detect available compute device."""
    if requested != "auto":
        return requested
    try:
        import torch
        if torch.cuda.is_available():
            gpu_name = torch.cuda.get_device_name(0)
            logger.info("CUDA available: %s", gpu_name)
            return "cuda"
        if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
            logger.info("MPS (Apple Silicon) available")
            return "mps"
    except ImportError:
        pass
    logger.info("Using CPU (no GPU detected)")
    return "cpu"


class RunPodClient:
    """Client for RunPod serverless GPU endpoints."""

    GRAPHQL_URL = "https://api.runpod.io/graphql"
    SERVERLESS_BASE = "https://api.runpod.ai/v2"

    def __init__(self, api_key: Optional[str] = None):
        self.api_key = api_key or os.environ.get("RUNPOD_API_KEY", "")
        if not self.api_key:
            logger.warning("RUNPOD_API_KEY not set; RunPod calls will fail")
        elif len(self.api_key) > 8:
            logger.info("RUNPOD_API_KEY loaded (ends with ...%s)", self.api_key[-4:])
        self.headers = {
            "Authorization": f"Bearer {self.api_key}",
            "Content-Type": "application/json",
        }

    def _sanitized_request(self, method: str, url: str, **kwargs: Any) -> Any:
        """Make a request with sanitized error logging (no API key leaks)."""
        import requests

        try:
            resp = getattr(requests, method)(url, headers=self.headers, **kwargs)
            resp.raise_for_status()
            return resp
        except requests.RequestException as exc:
            sanitized_msg = str(exc)
            if self.api_key and self.api_key in sanitized_msg:
                sanitized_msg = sanitized_msg.replace(self.api_key, "***REDACTED***")
            raise requests.RequestException(sanitized_msg) from None

    def submit_job(
        self,
        endpoint_id: str,
        payload: dict[str, Any],
        timeout: int = 60,
    ) -> str:
        """Submit a serverless job. Returns job ID."""
        url = f"{self.SERVERLESS_BASE}/{endpoint_id}/run"
        resp = self._sanitized_request("post", url, json={"input": payload}, timeout=timeout)
        data = resp.json()
        job_id = data.get("id", "")
        logger.info("RunPod job submitted: %s", job_id)
        return job_id

    def poll_job(
        self,
        endpoint_id: str,
        job_id: str,
        poll_interval: int = 10,
        max_wait: int = 3600,
    ) -> dict[str, Any]:
        """Poll a serverless job until completion."""
        url = f"{self.SERVERLESS_BASE}/{endpoint_id}/status/{job_id}"
        elapsed = 0
        while elapsed < max_wait:
            resp = self._sanitized_request("get", url, timeout=60)
            data = resp.json()
            status = data.get("status", "")
            if status == "COMPLETED":
                logger.info("RunPod job %s completed", job_id)
                return data.get("output", {})
            if status in ("FAILED", "CANCELLED"):
                raise RuntimeError(f"RunPod job {job_id} {status}: {data}")
            logger.debug("Job %s status: %s (elapsed %ds)", job_id, status, elapsed)
            time.sleep(poll_interval)
            elapsed += poll_interval
        raise TimeoutError(f"RunPod job {job_id} timed out after {max_wait}s")

    def run_job(
        self,
        endpoint_id: str,
        payload: dict[str, Any],
        poll_interval: int = 10,
        max_wait: int = 3600,
    ) -> dict[str, Any]:
        """Submit and poll a job in one call."""
        job_id = self.submit_job(endpoint_id, payload)
        return self.poll_job(endpoint_id, job_id, poll_interval, max_wait)
