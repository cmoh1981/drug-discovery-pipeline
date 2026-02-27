"""ZINC22 purchasable compound database API client.

ZINC22 (https://zinc22.docking.org) provides access to hundreds of millions
of purchasable drug-like molecules. It offers REST endpoints for similarity
and substructure search.

Note: ZINC22 API endpoints are rate-limited and may require an account for
large-scale queries. Free tier access supports moderate query volumes.
See https://zinc22.docking.org/doc/ for endpoint documentation.
"""

from __future__ import annotations

import logging
from typing import Any

from drugdiscovery.utils.web import get_with_retries, post_with_retries

logger = logging.getLogger(__name__)

ZINC_BASE = "https://zinc22.docking.org"

_HEADERS = {"Accept": "application/json"}


def search_zinc_by_similarity(
    smiles: str,
    threshold: float = 0.7,
    limit: int = 100,
) -> list[dict]:
    """Search ZINC22 for purchasable compounds by Tanimoto similarity.

    Endpoint: GET /substances/subsets/for-sale.json
              with params smiles, similarity, count

    Note: ZINC22's similarity threshold is expressed as a fraction (0.0–1.0).

    Args:
        smiles: Query SMILES string.
        threshold: Minimum Tanimoto similarity in [0, 1]. Values < 0.4 are
                   clamped to 0.4 to avoid excessively large result sets.
        limit: Maximum number of compounds to return.

    Returns:
        List of dicts with keys:
            id, smiles, name, source, vendor, purchasability
        Returns an empty list on any error.
    """
    try:
        threshold = max(0.4, min(1.0, float(threshold)))
        url = f"{ZINC_BASE}/substances/subsets/for-sale.json"
        params: dict[str, Any] = {
            "smiles": smiles,
            "similarity": threshold,
            "count": limit,
        }
        resp = get_with_retries(url, headers=_HEADERS, params=params)
        data = resp.json()
        results = _parse_substance_list(data, limit)
        logger.info(
            "ZINC22 search_by_similarity(threshold=%.2f): %d results",
            threshold,
            len(results),
        )
        return results

    except Exception as exc:  # noqa: BLE001
        logger.warning("ZINC22 search_by_similarity failed: %s", exc)
        return []


def search_zinc_by_substructure(smarts: str, limit: int = 100) -> list[dict]:
    """Search ZINC22 for purchasable compounds containing a SMARTS substructure.

    Endpoint: GET /substances/subsets/for-sale.json
              with params smarts, count

    Args:
        smarts: SMARTS pattern for substructure search.
        limit: Maximum number of compounds to return.

    Returns:
        List of dicts with keys:
            id, smiles, name, source, vendor, purchasability
        Returns an empty list on any error.
    """
    try:
        url = f"{ZINC_BASE}/substances/subsets/for-sale.json"
        params: dict[str, Any] = {"smarts": smarts, "count": limit}
        resp = get_with_retries(url, headers=_HEADERS, params=params)
        data = resp.json()
        results = _parse_substance_list(data, limit)
        logger.info(
            "ZINC22 search_by_substructure(%r): %d results", smarts, len(results)
        )
        return results

    except Exception as exc:  # noqa: BLE001
        logger.warning("ZINC22 search_by_substructure failed: %s", exc)
        return []


def search_zinc_by_name(name: str, limit: int = 50) -> list[dict]:
    """Search ZINC22 by compound name or catalog identifier.

    Endpoint: GET /substances.json?name={name}&count={limit}

    Args:
        name: Compound name or catalog ID (e.g. "ZINC000001234567").
        limit: Maximum number of results.

    Returns:
        Same schema as search_zinc_by_similarity().
        Returns an empty list on any error.
    """
    try:
        url = f"{ZINC_BASE}/substances.json"
        params: dict[str, Any] = {"name": name, "count": limit}
        resp = get_with_retries(url, headers=_HEADERS, params=params)
        data = resp.json()
        results = _parse_substance_list(data, limit)
        logger.info("ZINC22 search_by_name(%r): %d results", name, len(results))
        return results

    except Exception as exc:  # noqa: BLE001
        logger.warning("ZINC22 search_by_name failed: %s", exc)
        return []


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _parse_substance_list(data: dict | list, limit: int) -> list[dict]:
    """Normalise a ZINC22 API response to the pipeline standard list."""
    if isinstance(data, list):
        records = data
    else:
        records = (
            data.get("substances")
            or data.get("results")
            or data.get("data")
            or []
        )

    results: list[dict] = []
    for rec in records[:limit]:
        results.append(_substance_to_dict(rec))
    return results


def _substance_to_dict(rec: dict) -> dict:
    """Convert a raw ZINC22 substance record to the pipeline standard dict."""
    # ZINC IDs may appear under several field names across API versions
    zinc_id = (
        rec.get("zinc_id")
        or rec.get("ZINC_ID")
        or rec.get("id", "")
    )
    # Vendor / availability info
    vendor = ""
    purchasability = ""
    if "supplier" in rec:
        vendor = rec["supplier"]
    if "purchasability" in rec:
        purchasability = rec["purchasability"]
    elif "availability" in rec:
        purchasability = rec["availability"]

    return {
        "id": str(zinc_id),
        "smiles": rec.get("smiles") or rec.get("SMILES", ""),
        "name": rec.get("name") or rec.get("pref_name") or str(zinc_id),
        "source": "ZINC22",
        "vendor": vendor,
        "purchasability": purchasability,
    }
