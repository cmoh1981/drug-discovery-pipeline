"""Application settings via environment variables."""

from __future__ import annotations

from pathlib import Path

from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    """Service configuration loaded from .env or environment."""

    # Security
    SECRET_KEY: str = "change-me-in-production"
    ALGORITHM: str = "HS256"
    ACCESS_TOKEN_EXPIRE_MINUTES: int = 60 * 24  # 24 hours

    # Database
    DATABASE_URL: str = "sqlite:///./drugdiscovery_service.db"

    # Pipeline
    PIPELINE_OUTPUT_DIR: str = "job_results"
    MAX_CONCURRENT_JOBS: int = 2

    # RunPod
    RUNPOD_API_KEY: str = ""

    # CORS
    CORS_ORIGINS: str = "*"

    model_config = {"env_file": ".env", "env_file_encoding": "utf-8"}

    @property
    def output_path(self) -> Path:
        p = Path(self.PIPELINE_OUTPUT_DIR)
        p.mkdir(parents=True, exist_ok=True)
        return p


settings = Settings()
