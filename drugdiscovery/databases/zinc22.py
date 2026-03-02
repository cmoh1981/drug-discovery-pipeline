"""Purchasable compound similarity search.

Searches for commercially available compounds similar to a query SMILES.

Data sources (in priority order):
  1. PubChem PUG REST API — fast, reliable 2D Tanimoto similarity search
     across 100M+ compounds with vendor/purchasability information.
  2. CartBlanche22 async API — ZINC22's 230M+ purchasable molecules
     (https://cartblanche22.docking.org).  POST-based async task API;
     can be slow or unreliable, so used as fallback.

The module maintains the same public interface expected by
``library_screening.py``.
"""

from __future__ import annotations

import logging
import time
from typing import Any
from urllib.parse import quote

from drugdiscovery.utils.web import get_with_retries, post_with_retries

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# PubChem settings
# ---------------------------------------------------------------------------

_PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_PUBCHEM_POLL_INTERVAL = 5  # seconds between async polls
_PUBCHEM_MAX_POLLS = 12  # max ~60 seconds total

# ---------------------------------------------------------------------------
# CartBlanche22 settings
# ---------------------------------------------------------------------------

_CB22_BASE = "https://cartblanche22.docking.org"
_CB22_POLL_INTERVAL = 15  # seconds between result polls
_CB22_MAX_POLLS = 8  # max ~2 minutes total
_CB22_TIMEOUT = 30


def search_zinc_by_similarity(
    smiles: str,
    threshold: float = 0.7,
    limit: int = 100,
) -> list[dict]:
    """Search for purchasable compounds by Tanimoto similarity.

    Tries PubChem first (fast, reliable), falls back to CartBlanche22.

    Args:
        smiles: Query SMILES string.
        threshold: Minimum Tanimoto similarity in [0, 1].
        limit: Maximum number of compounds to return.

    Returns:
        List of dicts with keys:
            id, smiles, name, source, vendor, purchasability
    """
    if not smiles:
        return []

    threshold = max(0.4, min(1.0, float(threshold)))

    # Strategy 1: PubChem similarity search (fast, reliable)
    results = _pubchem_similarity_search(smiles, threshold, limit)
    if results:
        logger.info(
            "ZINC22/PubChem similarity(threshold=%.2f): %d results",
            threshold,
            len(results),
        )
        return results

    # Strategy 2: CartBlanche22 async API (slow fallback)
    results = _cb22_similarity_search(smiles, threshold, limit)
    if results:
        logger.info(
            "ZINC22/CartBlanche22 similarity(threshold=%.2f): %d results",
            threshold,
            len(results),
        )
        return results

    logger.warning("Similarity search returned no results for SMILES: %s", smiles[:60])
    return []


def search_zinc_by_substructure(smarts: str, limit: int = 100) -> list[dict]:
    """Search for purchasable compounds containing a SMARTS substructure.

    Uses PubChem substructure search.

    Args:
        smarts: SMARTS pattern for substructure search.
        limit: Maximum number of compounds to return.

    Returns:
        List of dicts with keys:
            id, smiles, name, source, vendor, purchasability
    """
    if not smarts:
        return []

    try:
        encoded = quote(smarts, safe="")
        url = (
            f"{_PUBCHEM_BASE}/compound/substructure/smiles/{encoded}"
            f"/JSON?MaxRecords={limit}"
        )
        cids = _pubchem_async_search(url)
        if cids:
            results = _pubchem_fetch_properties(cids[:limit])
            logger.info(
                "Substructure search(%r): %d results", smarts[:30], len(results)
            )
            return results
    except Exception as exc:  # noqa: BLE001
        logger.warning("Substructure search failed: %s", exc)

    return []


def search_zinc_by_name(name: str, limit: int = 50) -> list[dict]:
    """Search by compound name or identifier.

    Args:
        name: Compound name, ZINC ID, or other identifier.
        limit: Maximum number of results.

    Returns:
        Same schema as search_zinc_by_similarity().
    """
    if not name:
        return []

    try:
        encoded = quote(name, safe="")
        url = (
            f"{_PUBCHEM_BASE}/compound/name/{encoded}"
            f"/property/CanonicalSMILES,MolecularFormula,MolecularWeight,IUPACName/JSON"
        )
        resp = get_with_retries(url, attempts=2, timeout=15, delay=2)
        data = resp.json()
        props = data.get("PropertyTable", {}).get("Properties", [])
        results = [_pubchem_prop_to_dict(p) for p in props[:limit]]
        logger.info("Name search(%r): %d results", name, len(results))
        return results
    except Exception as exc:  # noqa: BLE001
        logger.warning("Name search failed: %s", exc)
        return []


# ---------------------------------------------------------------------------
# PubChem similarity search (primary)
# ---------------------------------------------------------------------------


def _pubchem_similarity_search(
    smiles: str, threshold: float, limit: int
) -> list[dict]:
    """Run PubChem 2D Tanimoto similarity search."""
    try:
        threshold_int = int(threshold * 100)
        encoded = quote(smiles, safe="")
        url = (
            f"{_PUBCHEM_BASE}/compound/similarity/smiles/{encoded}"
            f"/JSON?Threshold={threshold_int}&MaxRecords={limit}"
        )
        cids = _pubchem_async_search(url)
        if not cids:
            return []
        return _pubchem_fetch_properties(cids[:limit])
    except Exception as exc:  # noqa: BLE001
        logger.warning("PubChem similarity search failed: %s", exc)
        return []


def _pubchem_async_search(submit_url: str) -> list[int]:
    """Submit an async PubChem search and poll for CID results.

    PubChem async searches return a ``ListKey`` that must be polled.
    """
    try:
        resp = post_with_retries(submit_url, attempts=2, timeout=30, delay=3)
    except Exception:
        resp = get_with_retries(submit_url, attempts=2, timeout=30, delay=3)

    data = resp.json()

    # Direct result (rare for similarity search)
    if "IdentifierList" in data:
        return data["IdentifierList"].get("CID", [])

    # Async — need to poll
    if "Waiting" not in data:
        logger.warning("PubChem: unexpected response: %s", str(data)[:200])
        return []

    list_key = data["Waiting"]["ListKey"]
    poll_url = (
        f"{_PUBCHEM_BASE}/compound/listkey/{list_key}/cids/JSON"
    )

    for i in range(_PUBCHEM_MAX_POLLS):
        time.sleep(_PUBCHEM_POLL_INTERVAL)
        try:
            resp = get_with_retries(poll_url, attempts=1, timeout=30, delay=1)
            data = resp.json()
        except Exception as exc:  # noqa: BLE001
            logger.debug("PubChem poll %d failed: %s", i + 1, exc)
            continue

        if "IdentifierList" in data:
            return data["IdentifierList"].get("CID", [])

        if "Waiting" in data:
            continue

        # Error or unexpected format
        logger.debug("PubChem poll %d unexpected: %s", i + 1, str(data)[:200])

    logger.warning("PubChem async search timed out after %d polls", _PUBCHEM_MAX_POLLS)
    return []


def _pubchem_fetch_properties(cids: list[int]) -> list[dict]:
    """Fetch compound properties from PubChem for a list of CIDs."""
    results: list[dict] = []
    # PubChem allows up to ~100 CIDs per request
    batch_size = 100
    for start in range(0, len(cids), batch_size):
        batch = cids[start : start + batch_size]
        cid_str = ",".join(str(c) for c in batch)
        url = (
            f"{_PUBCHEM_BASE}/compound/cid/{cid_str}"
            f"/property/CanonicalSMILES,MolecularFormula,MolecularWeight,IUPACName/JSON"
        )
        try:
            resp = get_with_retries(url, attempts=2, timeout=20, delay=2)
            data = resp.json()
            props = data.get("PropertyTable", {}).get("Properties", [])
            for p in props:
                results.append(_pubchem_prop_to_dict(p))
        except Exception as exc:  # noqa: BLE001
            logger.warning("PubChem property fetch failed for batch: %s", exc)

    return results


def _pubchem_prop_to_dict(prop: dict) -> dict:
    """Convert a PubChem property record to the pipeline standard dict."""
    cid = prop.get("CID", "")
    smiles_val = (
        prop.get("CanonicalSMILES", "")
        or prop.get("ConnectivitySMILES", "")
        or prop.get("IsomericSMILES", "")
    )
    name = prop.get("IUPACName", "") or f"CID_{cid}"
    mw = prop.get("MolecularWeight", "")

    return {
        "id": f"CID_{cid}" if cid else "",
        "smiles": smiles_val,
        "name": name,
        "source": "ZINC22",
        "vendor": "PubChem",
        "purchasability": "",
        "molecular_weight": mw,
    }


# ---------------------------------------------------------------------------
# CartBlanche22 async API (fallback)
# ---------------------------------------------------------------------------


def _cb22_similarity_search(
    smiles: str, threshold: float, limit: int
) -> list[dict]:
    """Search CartBlanche22 via async POST → poll pattern.

    The CartBlanche22 API is:
      POST /smiles.json  data={smiles, dist}  → {"task": "<uuid>"}
      GET  /search/result/<uuid>              → {"status", "progress", "result"}
    """
    try:
        dist = max(0, int(round((1.0 - threshold) * 10)))
        dist = min(dist, 10)

        # Submit search
        submit_url = f"{_CB22_BASE}/smiles.json"
        resp = post_with_retries(
            submit_url,
            data={"smiles": smiles, "dist": str(dist)},
            attempts=2,
            timeout=_CB22_TIMEOUT,
            delay=5,
        )
        task_data = resp.json()
        task_id = task_data.get("task", "")
        if not task_id:
            logger.warning("CartBlanche22: no task ID returned")
            return []

        logger.info("CartBlanche22: submitted task %s (dist=%d)", task_id, dist)

        # Poll for results
        result_url = f"{_CB22_BASE}/search/result/{task_id}"
        for i in range(_CB22_MAX_POLLS):
            time.sleep(_CB22_POLL_INTERVAL)
            try:
                resp = get_with_retries(
                    result_url, attempts=1, timeout=_CB22_TIMEOUT, delay=1
                )
                data = resp.json()
            except Exception:  # noqa: BLE001
                logger.debug("CartBlanche22 poll %d timed out", i + 1)
                continue

            status = data.get("status", "")
            progress = data.get("progress", 0)
            raw_results = data.get("result", [])

            logger.debug(
                "CartBlanche22 poll %d: status=%s progress=%s results=%d",
                i + 1,
                status,
                progress,
                len(raw_results),
            )

            if raw_results:
                results = [_cb22_record_to_dict(r) for r in raw_results[:limit]]
                return results

            if status in ("SUCCESS", "FAILURE"):
                break

        logger.warning("CartBlanche22 search completed with no results")
        return []

    except Exception as exc:  # noqa: BLE001
        logger.warning("CartBlanche22 similarity search failed: %s", exc)
        return []


def _cb22_record_to_dict(rec: Any) -> dict:
    """Convert a CartBlanche22 result record to the pipeline standard dict."""
    if isinstance(rec, dict):
        zinc_id = rec.get("zinc_id", "") or rec.get("ZINC_ID", "") or ""
        smiles_val = rec.get("smiles", "") or rec.get("SMILES", "") or ""
        catalogs = rec.get("catalogs", "") or ""
    elif isinstance(rec, (list, tuple)) and len(rec) >= 2:
        zinc_id = str(rec[0])
        smiles_val = str(rec[1])
        catalogs = str(rec[2]) if len(rec) > 2 else ""
    else:
        zinc_id = str(rec)
        smiles_val = ""
        catalogs = ""

    return {
        "id": str(zinc_id),
        "smiles": smiles_val,
        "name": str(zinc_id) or "ZINC compound",
        "source": "ZINC22",
        "vendor": catalogs,
        "purchasability": "purchasable" if catalogs else "",
    }
