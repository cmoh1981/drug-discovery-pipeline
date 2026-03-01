"""SSE log streaming and progress parsing for pipeline jobs."""

from __future__ import annotations

import asyncio
import json
import re
from pathlib import Path

from service.config import settings

MODULE_STEPS: dict[str, tuple[str, int]] = {
    "[M1]": ("Target Preparation", 1),
    "[M2]": ("Library Screening", 2),
    "[M3]": ("De Novo Design", 3),
    "[M4]": ("Structure Prediction", 4),
    "[M4.5]": ("Molecular Docking", 5),
    "[M4.6]": ("Perturbation Biology", 6),
    "[M7]": ("ADMET Prediction", 7),
    "[M5]": ("Scoring", 8),
    "[M8]": ("Delivery System", 9),
    "[M9]": ("Final Report", 10),
}

_MODULE_PATTERN = re.compile(r"\[(M[\d.]+)\]")


def _find_log_file(job_id: str) -> Path | None:
    """Locate the pipeline.log for a job."""
    job_dir = settings.output_path / job_id
    if not job_dir.is_dir():
        return None
    logs = list(job_dir.rglob("pipeline.log"))
    return logs[0] if logs else None


def parse_module_tag(line: str) -> tuple[str, int] | None:
    """Extract module name and step number from a log line."""
    match = _MODULE_PATTERN.search(line)
    if match:
        tag = f"[{match.group(1)}]"
        if tag in MODULE_STEPS:
            return MODULE_STEPS[tag]
    return None


async def stream_logs(job_id: str):
    """Async generator that yields SSE events from pipeline.log.

    Tails the log file and parses module progress markers.
    """
    log_path = None
    # Wait up to 30s for the log file to appear
    for _ in range(60):
        log_path = _find_log_file(job_id)
        if log_path:
            break
        await asyncio.sleep(0.5)

    if not log_path:
        yield _sse_event({"type": "error", "message": "Log file not found"})
        return

    last_pos = 0
    idle_count = 0
    max_idle = 120  # stop after 60s of no new lines

    while idle_count < max_idle:
        try:
            with open(log_path, "r", encoding="utf-8", errors="replace") as f:
                f.seek(last_pos)
                new_lines = f.readlines()
                new_pos = f.tell()
        except FileNotFoundError:
            break

        if new_lines:
            idle_count = 0
            last_pos = new_pos
            for line in new_lines:
                line = line.rstrip()
                if not line:
                    continue
                event: dict = {"type": "log", "message": line}
                parsed = parse_module_tag(line)
                if parsed:
                    module_name, step = parsed
                    event["type"] = "progress"
                    event["module"] = module_name
                    event["step"] = step
                    event["total"] = 10
                yield _sse_event(event)

                # Check for pipeline completion
                if "PIPELINE COMPLETE" in line:
                    yield _sse_event({"type": "complete", "message": "Pipeline finished"})
                    return
        else:
            idle_count += 1
            await asyncio.sleep(0.5)

    yield _sse_event({"type": "timeout", "message": "Log stream timed out"})


def _sse_event(data: dict) -> str:
    """Format a dict as an SSE event string."""
    return f"data: {json.dumps(data)}\n\n"
