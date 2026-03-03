"""Subscription and payment management endpoints."""

from __future__ import annotations

import hashlib
import hmac
import logging
from datetime import datetime, timedelta, timezone

from fastapi import APIRouter, Depends, HTTPException, Request, status
from sqlalchemy.orm import Session

from service.config import settings
from service.database import get_db
from service.dependencies import get_current_user
from service.models import Job, Payment, Subscription, User
from service.portone import TIER_PRICES, get_portone_client
from service.schemas import (
    BillingKeyRequest,
    PaymentOut,
    QuotaInfo,
    SubscriptionCreateRequest,
    SubscriptionOut,
)

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/subscription", tags=["subscription"])

TIER_LIMITS = {
    "free": 3,
    "pro": 50,
    "enterprise": 999999,
}


def _utcnow() -> datetime:
    return datetime.now(timezone.utc)


def _get_active_subscription(db: Session, user_id: int) -> Subscription | None:
    return (
        db.query(Subscription)
        .filter(Subscription.user_id == user_id, Subscription.status.in_(("active", "past_due")))
        .first()
    )


def _count_total_jobs(db: Session, user_id: int) -> int:
    """Count total jobs for free-tier users (no subscription record)."""
    return db.query(Job).filter(Job.user_id == user_id).count()


# ── GET /api/subscription/config ─────────────────────────────────────

@router.get("/config")
def get_portone_config():
    """Return public PortOne keys for frontend SDK (no auth required)."""
    return {
        "store_id": settings.PORTONE_STORE_ID,
        "channel_key": settings.PORTONE_CHANNEL_KEY,
    }


# ── GET /api/subscription/ ──────────────────────────────────────────

@router.get("/", response_model=QuotaInfo)
def get_subscription(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    sub = _get_active_subscription(db, current_user.id)
    if sub and sub.tier != "free":
        tier = sub.tier
        runs_used = sub.runs_used_this_period
    else:
        tier = "free"
        runs_used = _count_total_jobs(db, current_user.id)

    runs_limit = TIER_LIMITS.get(tier, 3)
    return QuotaInfo(
        tier=tier,
        runs_used=runs_used,
        runs_limit=runs_limit,
        can_submit=runs_used < runs_limit,
        is_gpu_enabled=tier in ("pro", "enterprise"),
    )


# ── POST /api/subscription/register-billing-key ─────────────────────

@router.post("/register-billing-key")
def register_billing_key(
    body: BillingKeyRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    sub = _get_active_subscription(db, current_user.id)
    if not sub:
        sub = Subscription(user_id=current_user.id, tier="free", status="active")
        db.add(sub)
    sub.billing_key = body.billing_key
    db.commit()
    return {"status": "ok", "message": "Billing key registered"}


# ── POST /api/subscription/subscribe ────────────────────────────────

@router.post("/subscribe", response_model=SubscriptionOut)
async def subscribe(
    body: SubscriptionCreateRequest,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    sub = _get_active_subscription(db, current_user.id)
    if not sub or not sub.billing_key:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Register a billing key first",
        )

    if sub.tier == body.tier and sub.status == "active":
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Already subscribed to {body.tier}",
        )

    price = TIER_PRICES.get(body.tier)
    if price is None:
        raise HTTPException(status_code=400, detail="Invalid tier")

    # Charge first month via PortOne
    client = get_portone_client()
    try:
        result = await client.charge_billing_key(
            billing_key=sub.billing_key,
            amount=price,
            currency="KRW",
            order_name=f"Drug Discovery Pipeline - {body.tier.title()} Plan",
        )
    except Exception as e:
        logger.error("PortOne charge failed: %s", e)
        raise HTTPException(status_code=502, detail="Payment processing failed")

    # Verify the charged amount matches expected price (anti-fraud)
    payment_id = result.get("paymentId", "")
    try:
        verification = await client.verify_payment(payment_id)
        paid_amount = verification.get("amount", {}).get("total", 0)
        if paid_amount != price:
            logger.error("Amount mismatch: expected %d, got %d", price, paid_amount)
            await client.cancel_payment(payment_id, "Amount verification failed")
            raise HTTPException(status_code=400, detail="Payment amount verification failed")
    except HTTPException:
        raise
    except Exception as e:
        logger.error("Payment verification failed: %s", e)
        # Charge succeeded but verify failed -- proceed cautiously
        pass

    # Record payment
    now = _utcnow()
    payment = Payment(
        user_id=current_user.id,
        subscription_id=sub.id,
        portone_payment_id=payment_id,
        amount=price,
        currency="KRW",
        status="paid",
        tier=body.tier,
    )
    db.add(payment)

    # Update subscription
    sub.tier = body.tier
    sub.status = "active"
    sub.current_period_start = now
    sub.current_period_end = now + timedelta(days=30)
    sub.runs_used_this_period = 0
    sub.updated_at = now

    # Denormalize on user for fast access
    current_user.subscription_tier = body.tier
    db.commit()
    db.refresh(sub)
    return sub


# ── POST /api/subscription/cancel ───────────────────────────────────

@router.post("/cancel")
def cancel_subscription(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    sub = _get_active_subscription(db, current_user.id)
    if not sub or sub.tier == "free":
        raise HTTPException(status_code=400, detail="No active paid subscription")

    sub.status = "cancelled"
    sub.updated_at = _utcnow()
    # Keep access until period end -- tier stays until expiry
    db.commit()
    return {
        "status": "ok",
        "message": "Subscription cancelled. Access continues until period end.",
        "active_until": sub.current_period_end.isoformat() if sub.current_period_end else None,
    }


# ── GET /api/subscription/payments ──────────────────────────────────

@router.get("/payments", response_model=list[PaymentOut])
def list_payments(
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    payments = (
        db.query(Payment)
        .filter(Payment.user_id == current_user.id)
        .order_by(Payment.created_at.desc())
        .limit(50)
        .all()
    )
    return payments


# ── POST /api/subscription/webhook ──────────────────────────────────

@router.post("/webhook")
async def portone_webhook(request: Request, db: Session = Depends(get_db)):
    """Receive PortOne webhook notifications for payment status updates."""
    body_bytes = await request.body()
    # Verify webhook signature if secret is configured
    signature = request.headers.get("x-portone-signature")
    if signature and settings.PORTONE_API_SECRET:
        expected = hmac.new(
            settings.PORTONE_API_SECRET.encode(),
            body_bytes,
            hashlib.sha256,
        ).hexdigest()
        if not hmac.compare_digest(signature, expected):
            raise HTTPException(status_code=401, detail="Invalid webhook signature")

    import json
    try:
        data = json.loads(body_bytes)
    except json.JSONDecodeError:
        raise HTTPException(status_code=400, detail="Invalid JSON")

    event_type = data.get("type", "")
    payment_data = data.get("data", {})
    payment_id = payment_data.get("paymentId", "")

    if not payment_id:
        return {"status": "ignored"}

    payment = db.query(Payment).filter(Payment.portone_payment_id == payment_id).first()
    if not payment:
        logger.warning("Webhook for unknown payment: %s", payment_id)
        return {"status": "ignored"}

    if event_type == "Transaction.Paid":
        payment.status = "paid"
        if payment.subscription_id:
            sub = db.query(Subscription).filter(Subscription.id == payment.subscription_id).first()
            if sub:
                now = _utcnow()
                sub.current_period_end = now + timedelta(days=30)
                sub.runs_used_this_period = 0
                sub.status = "active"
                sub.updated_at = now
    elif event_type in ("Transaction.Failed", "Transaction.Cancelled"):
        payment.status = "failed"
        if payment.subscription_id:
            sub = db.query(Subscription).filter(Subscription.id == payment.subscription_id).first()
            if sub:
                sub.status = "past_due"
                sub.updated_at = _utcnow()

    db.commit()
    return {"status": "ok"}
