"""COCONUT (Collection of Open Natural Products) API client.

COCONUT is freely accessible at https://coconut.naturalproducts.net/
No authentication is required for the public search API.

API reference: https://coconut.naturalproducts.net/api
"""

from __future__ import annotations

import logging
from typing import Any

from drugdiscovery.utils.web import get_with_retries, post_with_retries

logger = logging.getLogger(__name__)

COCONUT_BASE = "https://coconut.naturalproducts.net/api"

_HEADERS = {"Accept": "application/json", "Content-Type": "application/json"}


def search_coconut(query: str, limit: int = 100) -> list[dict]:
    """Search COCONUT for natural products matching a text query.

    Sends a GET request to the simple search endpoint:
        GET /search/simple?query={query}&perPage={limit}

    Falls back to the POST /search endpoint if the simple endpoint is
    unavailable, as the COCONUT API has changed between versions.

    Args:
        query: Free-text query (compound name, structural class, taxonomy, etc.).
        limit: Maximum number of results to return.

    Returns:
        List of dicts with keys:
            id, smiles, name, source, molecular_formula, molecular_weight
        Returns an empty list on any error.
    """
    results = _search_simple(query, limit)
    if results is None:
        results = _search_post(query, limit)
    if results is None:
        logger.warning("COCONUT: both search strategies failed for query %r", query)
        return []

    logger.info("COCONUT search_coconut(%r): %d results", query, len(results))
    return results


def search_coconut_by_smiles(smiles: str, limit: int = 100) -> list[dict]:
    """Search COCONUT for natural products similar to a given SMILES.

    Uses the structure-based search endpoint:
        POST /search  with body {"type": "similarity", "smiles": smiles}

    Args:
        smiles: Query SMILES string.
        limit: Maximum number of results.

    Returns:
        Same schema as search_coconut().
        Returns an empty list on any error.
    """
    try:
        url = f"{COCONUT_BASE}/search"
        payload: dict[str, Any] = {
            "type": "similarity",
            "smiles": smiles,
            "perPage": limit,
        }
        resp = post_with_retries(url, json=payload, headers=_HEADERS)
        data = resp.json()
        results = _parse_response(data, limit)
        logger.info(
            "COCONUT search_coconut_by_smiles: %d results", len(results)
        )
        return results
    except Exception as exc:  # noqa: BLE001
        logger.warning("COCONUT search_coconut_by_smiles failed: %s", exc)
        return []


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _search_simple(query: str, limit: int) -> list[dict] | None:
    """GET /search/simple — returns None if the request fails."""
    try:
        url = f"{COCONUT_BASE}/search/simple"
        params: dict[str, Any] = {"query": query, "perPage": limit}
        resp = get_with_retries(url, headers=_HEADERS, params=params)
        data = resp.json()
        return _parse_response(data, limit)
    except Exception as exc:  # noqa: BLE001
        logger.debug("COCONUT simple search failed (will try POST): %s", exc)
        return None


def _search_post(query: str, limit: int) -> list[dict] | None:
    """POST /search with a text query — returns None if the request fails."""
    try:
        url = f"{COCONUT_BASE}/search"
        payload: dict[str, Any] = {
            "type": "text",
            "query": query,
            "perPage": limit,
        }
        resp = post_with_retries(url, json=payload, headers=_HEADERS)
        data = resp.json()
        return _parse_response(data, limit)
    except Exception as exc:  # noqa: BLE001
        logger.debug("COCONUT POST search failed: %s", exc)
        return None


def _parse_response(data: dict | list, limit: int) -> list[dict]:
    """Normalise a COCONUT API response to the pipeline standard list."""
    # The API may return {"data": [...]} or a bare list depending on version
    if isinstance(data, list):
        records = data
    else:
        records = (
            data.get("data")
            or data.get("results")
            or data.get("molecules")
            or []
        )

    results: list[dict] = []
    for rec in records[:limit]:
        results.append(_record_to_dict(rec))
    return results


def _record_to_dict(rec: dict) -> dict:
    """Convert a raw COCONUT record to the pipeline standard dict."""
    # COCONUT uses 'coconut_id' or 'identifier' depending on API version
    identifier = (
        rec.get("coconut_id")
        or rec.get("identifier")
        or rec.get("id", "")
    )
    return {
        "id": str(identifier),
        "smiles": rec.get("smiles") or rec.get("canonical_smiles", ""),
        "name": rec.get("name") or rec.get("iupac_name", ""),
        "source": "COCONUT",
        "molecular_formula": rec.get("molecular_formula", ""),
        "molecular_weight": rec.get("molecular_weight") or rec.get("exact_molecular_weight"),
    }
