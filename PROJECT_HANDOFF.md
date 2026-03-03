# Drug Discovery Pipeline — Project Handoff Guide

## Quick Start (New Laptop)

### 1. Prerequisites
```bash
# Python 3.10+
python --version

# Install dependencies
pip install -e ".[web]"

# Or install from requirements.txt
pip install -r requirements.txt
pip install fastapi uvicorn sqlalchemy python-jose[cryptography] bcrypt python-multipart pydantic-settings email-validator httpx
```

### 2. Environment Setup
```bash
# Copy env template
cp .env.example .env

# Edit .env — set at minimum:
SECRET_KEY=your-random-secret-key-here
```

### 3. Run Locally
```bash
python run_server.py
# → http://localhost:8000
# → Swagger docs: http://localhost:8000/docs
```

### 4. Test Account
- Email: `test@test.com`
- Password: `test`

---

## Architecture

### Project Structure
```
drug-discovery-pipeline/
├── drugdiscovery/           # Core pipeline (10 modules: M1-M9 + M4.5)
│   ├── modules/             # Pipeline modules (target_prep, docking, admet, etc.)
│   ├── config.py            # PipelineConfig builder
│   ├── types.py             # Candidate, TargetProfile dataclasses
│   └── io.py                # CSV/JSON I/O utilities
├── service/                 # FastAPI web service (SaaS layer)
│   ├── app.py               # Application factory + CORS + lifespan
│   ├── config.py            # Settings via pydantic-settings (.env)
│   ├── database.py          # SQLAlchemy engine (SQLite WAL)
│   ├── models.py            # ORM: User, Job, Subscription, Payment
│   ├── schemas.py           # Pydantic V2 request/response models
│   ├── auth.py              # JWT + bcrypt password hashing
│   ├── portone.py           # PortOne V2 payment API client (httpx)
│   ├── scheduler.py         # Monthly billing scheduler (background)
│   ├── worker.py            # ThreadPoolExecutor for pipeline jobs
│   ├── routers/
│   │   ├── auth_router.py        # /api/auth/* (register, login, profile)
│   │   ├── jobs_router.py        # /api/jobs/* (submit, list, progress, cancel)
│   │   ├── results_router.py     # /api/results/* (candidates, dashboard, pubchem)
│   │   ├── compare_router.py     # /api/compare/* (cross-run comparison)
│   │   ├── subscription_router.py # /api/subscription/* (plans, payments, webhook)
│   │   └── health_router.py      # /api/health
│   └── static/              # Frontend SPA (vanilla JS, no bundler)
│       ├── index.html        # SPA shell with hash routing
│       ├── css/
│       │   ├── design-tokens.css  # Color palette, typography
│       │   ├── layout.css         # Sidebar, top bar, responsive
│       │   ├── components.css     # Buttons, cards, forms, badges
│       │   ├── pages.css          # Page-specific styles
│       │   └── pricing.css        # Pricing page styles
│       └── js/
│           ├── api.js             # Auth-aware fetch wrapper + SSE
│           ├── router.js          # Hash-based router with :param support
│           ├── app.js             # Route registration, theme, auth guard
│           └── pages/
│               ├── landing.js      # Home page (#/)
│               ├── login.js        # Sign in (#/login)
│               ├── register.js     # Create account (#/register)
│               ├── jobs-gallery.js  # Job list (#/jobs)
│               ├── job-submit.js    # New job form (#/jobs/new)
│               ├── job-detail.js    # Job detail (#/jobs/:id)
│               ├── job-dashboard.js # Dashboard iframe (#/jobs/:id/dashboard)
│               ├── compare.js       # Compare runs (#/compare)
│               └── pricing.js       # Pricing & subscription (#/pricing)
├── tests/                   # Pipeline unit tests (~197)
├── tests_service/           # Service API tests (~28)
├── configs/                 # YAML config templates
├── Dockerfile               # Docker build for Railway
├── pyproject.toml           # Python package config
├── requirements.txt         # Pip dependencies
├── run_server.py            # Uvicorn launcher
└── railway.toml             # Railway deployment config
```

### Pipeline Flow (10 Modules)
```
M1 Target Prep → M2 Library Screen → M3 De Novo Design →
M4 Structure Prediction → M4.5 Molecular Docking →
M4.6 Perturbation Biology → M5 Scoring →
M7 ADMET → M8 Delivery → M9 Report
```

### Subscription Tiers (PortOne V2 / Korean Won)
| Tier | Price | Runs | GPU |
|------|-------|------|-----|
| Free | ₩0 | 3 total | No |
| Pro | ₩39,000/월 | 50/month | Yes |
| Enterprise | ₩129,000/월 | Unlimited | Yes |

---

## Deployment (Railway)

### Current Setup
- **Platform**: Railway (https://railway.com)
- **URL**: https://drug-discovery-pipeline-production.up.railway.app
- **GitHub**: https://github.com/cmoh1981/drug-discovery-pipeline
- **Build**: Dockerfile-based

### Deploy Commands
```bash
# Push to GitHub (auto-deploys if configured)
git push origin master

# Manual deploy via Railway CLI
railway up

# Check logs
railway logs
```

### Railway Environment Variables
```
SECRET_KEY=<your-jwt-secret>
PORTONE_API_SECRET=<portone-v2-secret>
PORTONE_STORE_ID=<portone-store-id>
PORTONE_CHANNEL_KEY=<portone-channel-key>
PORT=8000
```

---

## Key Technologies
- **Backend**: FastAPI, SQLAlchemy (SQLite WAL), JWT auth (python-jose), bcrypt
- **Frontend**: Vanilla JS SPA, hash routing, no bundler
- **Payments**: PortOne V2 REST API (httpx), billing key flow for subscriptions
- **Pipeline**: RDKit, Biopython, AutoDock Vina, PepMLM (optional GPU)
- **Deployment**: Docker on Railway

## API Endpoints Summary
| Group | Endpoints |
|-------|-----------|
| Auth | POST /api/auth/register, /login, GET /me |
| Jobs | POST /api/jobs/, GET list, GET /:id, GET /:id/progress, DELETE /:id |
| Results | GET /api/results/:id, /candidates, /dashboard, /pubchem |
| Compare | POST /api/compare |
| Subscription | GET /api/subscription/, POST /subscribe, /cancel, /webhook |
| Health | GET /api/health |

## Tests
```bash
# All tests
python -m pytest tests/ tests_service/ -v

# Service only
python -m pytest tests_service/ -v
```

## Recent Changes (March 2025)
1. PortOne subscription system (billing key flow, 3 tiers, quota enforcement)
2. Pricing page accessible without login (service info, FAQ, pipeline overview)
3. UI refresh: bright teal palette, light sidebar, clickable brand
4. Prices switched to Korean Won (₩)
5. Dashboard iframe URL fix (Windows backslash paths)
6. PortOne SDK async loading fix
