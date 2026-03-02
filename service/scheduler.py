"""Background scheduler for recurring subscription billing."""

from __future__ import annotations

import asyncio
import logging
from datetime import datetime, timedelta, timezone

from service.database import SessionLocal
from service.models import Payment, Subscription
from service.portone import TIER_PRICES, get_portone_client

logger = logging.getLogger(__name__)

MAX_RETRY_DAYS = 7
CHECK_INTERVAL_SECONDS = 3600  # 1 hour


def _utcnow() -> datetime:
    return datetime.now(timezone.utc)


async def _process_renewal(sub: Subscription) -> None:
    """Attempt to charge a single subscription renewal."""
    db = SessionLocal()
    try:
        # Re-fetch inside fresh session
        sub = db.query(Subscription).filter(Subscription.id == sub.id).first()
        if not sub or not sub.billing_key:
            return

        price = TIER_PRICES.get(sub.tier)
        if not price:
            return

        client = get_portone_client()
        try:
            result = await client.charge_billing_key(
                billing_key=sub.billing_key,
                amount=price,
                currency="USD",
                order_name=f"Drug Discovery Pipeline - {sub.tier.title()} Plan (Renewal)",
            )
            payment_id = result.get("paymentId", "")

            payment = Payment(
                user_id=sub.user_id,
                subscription_id=sub.id,
                portone_payment_id=payment_id,
                amount=price,
                currency="USD",
                status="paid",
                tier=sub.tier,
            )
            db.add(payment)

            now = _utcnow()
            sub.current_period_start = now
            sub.current_period_end = now + timedelta(days=30)
            sub.runs_used_this_period = 0
            sub.status = "active"
            sub.updated_at = now
            db.commit()
            logger.info("Renewed subscription %s for user %s", sub.id, sub.user_id)

        except Exception as e:
            logger.error("Renewal charge failed for sub %s: %s", sub.id, e)
            now = _utcnow()
            sub.status = "past_due"
            sub.updated_at = now

            # Expire after MAX_RETRY_DAYS of failed attempts
            if sub.current_period_end and (now - sub.current_period_end).days > MAX_RETRY_DAYS:
                sub.status = "expired"
                from service.models import User
                user = db.query(User).filter(User.id == sub.user_id).first()
                if user:
                    user.subscription_tier = "free"
                logger.warning("Subscription %s expired after retry period", sub.id)

            db.commit()
    finally:
        db.close()


async def billing_check_loop() -> None:
    """Periodically check for subscriptions needing renewal."""
    logger.info("Billing scheduler started (interval=%ds)", CHECK_INTERVAL_SECONDS)
    while True:
        try:
            db = SessionLocal()
            now = _utcnow()
            expired_subs = (
                db.query(Subscription)
                .filter(
                    Subscription.status.in_(("active", "past_due")),
                    Subscription.tier.in_(("pro", "enterprise")),
                    Subscription.billing_key.isnot(None),
                    Subscription.current_period_end < now,
                )
                .all()
            )
            db.close()

            for sub in expired_subs:
                await _process_renewal(sub)

            if expired_subs:
                logger.info("Processed %d subscription renewals", len(expired_subs))

        except Exception as e:
            logger.error("Billing check error: %s", e)

        await asyncio.sleep(CHECK_INTERVAL_SECONDS)
