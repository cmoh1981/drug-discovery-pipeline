"""PortOne V2 REST API client using httpx."""

from __future__ import annotations

import logging
import uuid

import httpx

logger = logging.getLogger(__name__)

TIER_PRICES = {
    "pro": 2900,        # $29.00 in cents
    "enterprise": 9900,  # $99.00 in cents
}


class PortOneClient:
    """Async wrapper around the PortOne V2 payment API."""

    BASE_URL = "https://api.portone.io"

    def __init__(self, api_secret: str, store_id: str):
        self.api_secret = api_secret
        self.store_id = store_id
        self.headers = {
            "Authorization": f"PortOne {api_secret}",
            "Content-Type": "application/json",
        }

    async def verify_payment(self, payment_id: str) -> dict:
        """GET /payments/{payment_id} -- verify payment status & amount."""
        async with httpx.AsyncClient() as client:
            resp = await client.get(
                f"{self.BASE_URL}/payments/{payment_id}",
                headers=self.headers,
            )
            resp.raise_for_status()
            return resp.json()

    async def charge_billing_key(
        self,
        billing_key: str,
        amount: int,
        currency: str = "USD",
        order_name: str = "Drug Discovery Pipeline Subscription",
    ) -> dict:
        """POST /payments/{payment_id}/billing-key -- charge a saved card."""
        payment_id = f"pay_{uuid.uuid4().hex[:20]}"
        payload = {
            "storeId": self.store_id,
            "billingKey": billing_key,
            "orderName": order_name,
            "amount": {"total": amount, "currency": currency},
        }
        async with httpx.AsyncClient() as client:
            resp = await client.post(
                f"{self.BASE_URL}/payments/{payment_id}/billing-key",
                headers=self.headers,
                json=payload,
            )
            resp.raise_for_status()
            data = resp.json()
            data["paymentId"] = payment_id
            return data

    async def get_billing_key_info(self, billing_key: str) -> dict:
        """GET /billing-keys/{billing_key} -- get saved card info."""
        async with httpx.AsyncClient() as client:
            resp = await client.get(
                f"{self.BASE_URL}/billing-keys/{billing_key}",
                headers=self.headers,
            )
            resp.raise_for_status()
            return resp.json()

    async def cancel_payment(self, payment_id: str, reason: str) -> dict:
        """POST /payments/{payment_id}/cancel -- refund a payment."""
        async with httpx.AsyncClient() as client:
            resp = await client.post(
                f"{self.BASE_URL}/payments/{payment_id}/cancel",
                headers=self.headers,
                json={"reason": reason},
            )
            resp.raise_for_status()
            return resp.json()


def get_portone_client() -> PortOneClient:
    """Create a PortOneClient from application settings."""
    from service.config import settings

    return PortOneClient(
        api_secret=settings.PORTONE_API_SECRET,
        store_id=settings.PORTONE_STORE_ID,
    )
