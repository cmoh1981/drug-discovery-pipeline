"""FastAPI application factory."""

from __future__ import annotations

import asyncio
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles

from service.config import settings
from service.database import create_tables
from service.worker import job_manager


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Startup and shutdown events."""
    create_tables()
    settings.output_path.mkdir(parents=True, exist_ok=True)

    # Start billing scheduler if PortOne is configured
    scheduler_task = None
    if settings.PORTONE_API_SECRET:
        from service.scheduler import billing_check_loop
        scheduler_task = asyncio.create_task(billing_check_loop())

    yield

    if scheduler_task:
        scheduler_task.cancel()
    job_manager.shutdown(wait=False)


def create_app() -> FastAPI:
    app = FastAPI(
        title="Drug Discovery Pipeline API",
        description="SaaS platform for automated drug discovery",
        version="1.0.0",
        lifespan=lifespan,
    )

    # CORS
    origins = [o.strip() for o in settings.CORS_ORIGINS.split(",") if o.strip()]
    app.add_middleware(
        CORSMiddleware,
        allow_origins=origins,
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    # Routers
    from service.routers.auth_router import router as auth_router
    from service.routers.jobs_router import router as jobs_router
    from service.routers.health_router import router as health_router

    app.include_router(auth_router)
    app.include_router(jobs_router)
    app.include_router(health_router)

    # Phase 2+ routers (imported conditionally to allow incremental deployment)
    try:
        from service.routers.results_router import router as results_router
        from service.routers.results_router import pubchem_router
        app.include_router(results_router)
        app.include_router(pubchem_router)
    except ImportError:
        pass

    try:
        from service.routers.compare_router import router as compare_router
        app.include_router(compare_router)
    except ImportError:
        pass

    try:
        from service.routers.subscription_router import router as subscription_router
        app.include_router(subscription_router)
    except Exception as exc:
        import logging
        logging.getLogger(__name__).warning("Subscription router not loaded: %s", exc)

    # Mount static files for pipeline reports
    import os
    reports_dir = settings.PIPELINE_OUTPUT_DIR
    if os.path.isdir(reports_dir):
        app.mount("/reports", StaticFiles(directory=reports_dir), name="reports")

    # Mount web UI static assets
    static_dir = Path(__file__).parent / "static"
    if static_dir.is_dir():
        app.mount("/static", StaticFiles(directory=str(static_dir)), name="static_assets")

    # Serve SPA root — hash routing means the browser always requests /
    index_html = static_dir / "index.html"

    @app.get("/", response_class=HTMLResponse)
    def serve_spa():
        return index_html.read_text(encoding="utf-8")

    return app
