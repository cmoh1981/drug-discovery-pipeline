# Update 2026-03-02 — PortOne Subscription & Payment System

## Summary

Added a complete PortOne V2 payment gateway integration with tiered subscription plans to monetize GPU access and enforce pipeline run quotas.

## Subscription Tiers

| Tier | Price | Runs | GPU (RunPod) |
|------|-------|------|--------------|
| Free | $0 | 3 total | No |
| Pro | $29/mo | 50/month | Yes |
| Enterprise | $99/mo | Unlimited | Yes |

## New Files (5)

| File | Description |
|------|-------------|
| `service/portone.py` | PortOne V2 API client (httpx) — verify, charge, cancel payments |
| `service/routers/subscription_router.py` | 7 endpoints: quota info, billing key registration, subscribe, cancel, payment history, webhook, config |
| `service/scheduler.py` | Background billing loop — auto-renews subscriptions, handles failures, expires after 7 days past-due |
| `service/static/js/pages/pricing.js` | Pricing page with 3-column plan cards, PortOne SDK popup, cancel flow |
| `service/static/css/pricing.css` | Responsive pricing grid, featured card highlight, payment history table |

## Modified Files (9)

| File | Changes |
|------|---------|
| `service/models.py` | Added `Subscription` and `Payment` ORM models; added `subscription_tier` column + relationships to `User` |
| `service/schemas.py` | Added `QuotaInfo`, `SubscriptionOut`, `SubscriptionCreateRequest`, `BillingKeyRequest`, `PaymentOut` Pydantic schemas |
| `service/config.py` | Added `PORTONE_API_SECRET`, `PORTONE_STORE_ID`, `PORTONE_CHANNEL_KEY` settings |
| `service/app.py` | Registered subscription router; starts billing scheduler on startup if PortOne is configured |
| `service/routers/jobs_router.py` | Added quota check before job creation (402 if exceeded); blocks GPU on free tier; increments `runs_used_this_period` for paid users |
| `service/database.py` | Added `_run_migrations()` to auto-add `subscription_tier` column to existing `users` table on startup |
| `service/static/index.html` | Added PortOne SDK script, pricing.css link, pricing.js script, "Pricing" nav link in sidebar |
| `service/static/js/app.js` | Registered `#/pricing` route; loads PortOne config on boot; added pricing to page titles |
| `service/static/js/pages/job-submit.js` | Added quota display above submit button; disables submit when quota exceeded; blocks GPU checkbox on free tier |

## API Endpoints Added

| Method | Path | Auth | Description |
|--------|------|------|-------------|
| GET | `/api/subscription/config` | No | Public PortOne keys for frontend SDK |
| GET | `/api/subscription/` | Yes | Current user quota info (tier, runs used/limit, can_submit, gpu_enabled) |
| POST | `/api/subscription/register-billing-key` | Yes | Store billing key from PortOne SDK popup |
| POST | `/api/subscription/subscribe` | Yes | Subscribe to Pro/Enterprise (charges first month) |
| POST | `/api/subscription/cancel` | Yes | Cancel subscription (access until period end) |
| GET | `/api/subscription/payments` | Yes | Payment history (last 50) |
| POST | `/api/subscription/webhook` | No | PortOne webhook receiver (HMAC-SHA256 verified) |

## Security

- Server-side payment amount verification (anti-fraud)
- HMAC-SHA256 webhook signature validation
- Billing key stored server-side only
- Quota enforcement is server-side (cannot be bypassed from frontend)

## Environment Variables (add to Railway)

```
PORTONE_API_SECRET=    # PortOne V2 API Secret
PORTONE_STORE_ID=      # PortOne Store ID
PORTONE_CHANNEL_KEY=   # Payment channel key (frontend SDK)
```

## Verification

- All 224 existing tests pass (zero regressions)
- E2E flow confirmed: register → free tier (3 runs) → 4th job blocked (402) → GPU blocked on free tier
- Database migration handles both fresh installs and existing SQLite databases
