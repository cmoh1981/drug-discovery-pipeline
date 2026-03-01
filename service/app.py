"""FastAPI application factory."""

from __future__ import annotations

from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles

from service.config import settings
from service.database import create_tables
from service.worker import job_manager


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Startup and shutdown events."""
    create_tables()
    settings.output_path.mkdir(parents=True, exist_ok=True)
    yield
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
        app.include_router(results_router)
    except ImportError:
        pass

    try:
        from service.routers.compare_router import router as compare_router
        app.include_router(compare_router)
    except ImportError:
        pass

    # Mount static files for pipeline reports
    import os
    reports_dir = settings.PIPELINE_OUTPUT_DIR
    if os.path.isdir(reports_dir):
        app.mount("/reports", StaticFiles(directory=reports_dir), name="reports")

    return app
