"""Pydantic V2 request/response schemas."""

from __future__ import annotations

from datetime import datetime
from typing import Any

from pydantic import BaseModel, EmailStr, Field


# ── Auth ──────────────────────────────────────────────────────────────

class UserCreate(BaseModel):
    email: EmailStr
    password: str = Field(min_length=8)
    full_name: str = ""
    organization: str = ""


class UserResponse(BaseModel):
    id: int
    email: str
    full_name: str
    organization: str
    is_active: bool
    created_at: datetime

    model_config = {"from_attributes": True}


class Token(BaseModel):
    access_token: str
    token_type: str = "bearer"


class LoginRequest(BaseModel):
    email: EmailStr
    password: str


# ── Jobs ──────────────────────────────────────────────────────────────

class JobCreate(BaseModel):
    modality: str = Field(pattern=r"^(small_molecule|peptide)$")
    mode: str = Field(pattern=r"^(agonist|antagonist)$")
    target: str = Field(min_length=1, max_length=255)
    tissue: str = ""
    top_n: int = Field(default=20, ge=1, le=500)
    use_runpod: bool = False
    config_overrides: dict[str, Any] | None = None


class JobResponse(BaseModel):
    id: str
    user_id: int
    status: str
    modality: str
    mode: str
    target: str
    tissue: str
    top_n: int
    use_runpod: bool
    config_overrides: str | None = None
    result_dir: str | None = None
    error_message: str | None = None
    total_candidates: int | None = None
    top_score: float | None = None
    current_module: str | None = None
    progress_step: int = 0
    progress_total: int = 10
    started_at: datetime | None = None
    completed_at: datetime | None = None
    created_at: datetime

    model_config = {"from_attributes": True}


class JobListResponse(BaseModel):
    jobs: list[JobResponse]
    total: int
    page: int
    page_size: int


# ── Results ───────────────────────────────────────────────────────────

class CandidateResponse(BaseModel):
    candidate_id: str
    candidate_type: str = ""
    modality: str = ""
    source: str = ""
    smiles: str = ""
    iupac_name: str = ""
    sequence: str = ""
    molecular_weight: float = 0.0
    binding_score: float = 0.0
    selectivity_score: float = 0.0
    drug_likeness: float = 0.0
    admet_score: float = 0.0
    admet_flags: int = 0
    delivery_system: str = ""
    perturbation_score: float = 0.0
    composite_score: float = 0.0
    rank: int = 0


class CandidateListResponse(BaseModel):
    candidates: list[CandidateResponse]
    total: int
    page: int
    page_size: int


class JobSummary(BaseModel):
    job_id: str
    status: str
    target: str
    modality: str
    mode: str
    total_candidates: int | None = None
    top_score: float | None = None
    files: list[str] = []
    started_at: datetime | None = None
    completed_at: datetime | None = None


class ProgressResponse(BaseModel):
    step: int
    total: int
    module: str | None = None
    status: str


# ── Comparison ────────────────────────────────────────────────────────

class CompareRequest(BaseModel):
    job_ids: list[str] = Field(min_length=2, max_length=10)
    top_n: int = Field(default=20, ge=1, le=100)


class ComparedCandidate(BaseModel):
    identifier: str
    scores: dict[str, float | None]


class CompareResponse(BaseModel):
    job_ids: list[str]
    shared_candidates: list[ComparedCandidate]
    unique_per_job: dict[str, int]
    total_compared: int


# ── Templates ─────────────────────────────────────────────────────────

class TemplateCreate(BaseModel):
    name: str = Field(min_length=1, max_length=255)
    description: str = ""
    config: dict[str, Any]


class TemplateResponse(BaseModel):
    id: int
    name: str
    description: str
    config_json: str
    created_at: datetime

    model_config = {"from_attributes": True}


# ── Health ────────────────────────────────────────────────────────────

class HealthResponse(BaseModel):
    status: str = "ok"
    database: str = "connected"
    active_jobs: int = 0
