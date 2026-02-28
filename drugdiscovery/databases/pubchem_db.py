"""PubChem Bioassay database client for target-based compound discovery.

PubChem PUG REST API provides access to bioassay data, allowing searches
for compounds with measured bioactivity against specific gene targets.

API reference: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest
Rate limit: 5 requests/second (we use a 0.25s delay between requests).
"""

from __future__ import annotations

import logging
import time
from typing import Any

from drugdiscovery.utils.web import get_with_retries

logger = logging.getLogger(__name__)

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

_HEADERS = {"Accept": "application/json"}

# PubChem enforces 5 requests/second; 0.25s delay keeps us within budget.
_RATE_DELAY = 0.25


def search_by_target(gene_name: str, limit: int = 200) -> list[dict]:
    """Search PubChem for compounds with bioactivity against a gene target.

    Workflow:
      1. Find bioassay IDs (AIDs) associated with the gene symbol.
      2. For each assay, retrieve CIDs of compounds flagged as 'active'.
      3. Batch-fetch compound properties (SMILES, MW, IUPAC name).

    Args:
        gene_name: HGNC gene symbol (e.g. ``"EGFR"``, ``"BRAF"``).
        limit: Maximum number of compounds to return across all assays.

    Returns:
        List of dicts with keys:
            id, smiles, name, source, molecular_weight,
            activity_type, activity_value
        Returns an empty list on any error.
    """
    try:
        aids = _get_assay_ids(gene_name, max_assays=10)
        if not aids:
            logger.warning(
                "PubChem: no bioassays found for gene %r", gene_name
            )
            return []

        logger.info("PubChem: found %d assays for gene %r", len(aids), gene_name)

        # Collect active compound CIDs across all assays
        all_cids: list[int] = []
        for aid in aids:
            if len(all_cids) >= limit:
                break
            time.sleep(_RATE_DELAY)
            cids = _get_active_cids(aid, max_cids=50)
            all_cids.extend(cids)

        if not all_cids:
            logger.info("PubChem: no active compounds for gene %r", gene_name)
            return []

        # Deduplicate and apply limit
        unique_cids = list(dict.fromkeys(all_cids))[:limit]
        logger.info(
            "PubChem: %d unique active CIDs (from %d total) for gene %r",
            len(unique_cids),
            len(all_cids),
            gene_name,
        )

        # Fetch compound properties in batches
        results = _get_compound_properties(unique_cids)
        logger.info(
            "PubChem search_by_target(%r): %d results", gene_name, len(results)
        )
        return results

    except Exception as exc:  # noqa: BLE001
        logger.warning("PubChem search_by_target failed: %s", exc)
        return []


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _get_assay_ids(gene_name: str, max_assays: int = 10) -> list[int]:
    """Retrieve bioassay IDs (AIDs) linked to a gene symbol.

    Endpoint:
        GET /assay/target/genesymbol/{gene}/aids/JSON

    Args:
        gene_name: HGNC gene symbol.
        max_assays: Maximum number of assay IDs to return.

    Returns:
        List of integer assay IDs, or empty list on failure.
    """
    url = f"{PUBCHEM_BASE}/assay/target/genesymbol/{gene_name}/aids/JSON"
    try:
        resp = get_with_retries(url, headers=_HEADERS, timeout=30, attempts=2, delay=3)
        data: dict[str, Any] = resp.json()
        aids = (
            data
            .get("InformationList", {})
            .get("Information", [{}])[0]
            .get("AID", [])
        )
        return [int(a) for a in aids[:max_assays]]
    except Exception as exc:  # noqa: BLE001
        logger.debug("PubChem _get_assay_ids(%r) failed: %s", gene_name, exc)
        return []


def _get_active_cids(aid: int, max_cids: int = 50) -> list[int]:
    """Retrieve CIDs of compounds flagged as 'active' in a bioassay.

    Endpoint:
        GET /assay/aid/{aid}/cids/JSON?cids_type=active

    Args:
        aid: PubChem bioassay ID.
        max_cids: Maximum number of CIDs to return per assay.

    Returns:
        List of integer compound IDs, or empty list on failure.
    """
    url = f"{PUBCHEM_BASE}/assay/aid/{aid}/cids/JSON"
    params: dict[str, str] = {"cids_type": "active"}
    try:
        resp = get_with_retries(
            url, headers=_HEADERS, params=params, timeout=30, attempts=2, delay=3
        )
        data: dict[str, Any] = resp.json()
        cids = (
            data
            .get("InformationList", {})
            .get("Information", [{}])[0]
            .get("CID", [])
        )
        return [int(c) for c in cids[:max_cids]]
    except Exception as exc:  # noqa: BLE001
        logger.debug("PubChem _get_active_cids(aid=%d) failed: %s", aid, exc)
        return []


def _get_compound_properties(cids: list[int]) -> list[dict]:
    """Batch-fetch compound properties from PubChem.

    Endpoint:
        GET /compound/cid/{cid_csv}/property/CanonicalSMILES,MolecularWeight,IUPACName/JSON

    CIDs are batched in groups of 100 to respect URL length limits.

    Args:
        cids: List of PubChem compound IDs.

    Returns:
        List of pipeline-standard dicts with compound data.
    """
    results: list[dict] = []
    batch_size = 100

    for i in range(0, len(cids), batch_size):
        batch = cids[i : i + batch_size]
        cid_csv = ",".join(str(c) for c in batch)
        url = (
            f"{PUBCHEM_BASE}/compound/cid/{cid_csv}"
            f"/property/CanonicalSMILES,MolecularWeight,IUPACName/JSON"
        )

        time.sleep(_RATE_DELAY)

        try:
            resp = get_with_retries(
                url, headers=_HEADERS, timeout=30, attempts=2, delay=3
            )
            data: dict[str, Any] = resp.json()
            properties = data.get("PropertyTable", {}).get("Properties", [])
            for prop in properties:
                results.append(_property_to_dict(prop))
        except Exception as exc:  # noqa: BLE001
            logger.debug(
                "PubChem _get_compound_properties batch %d-%d failed: %s",
                i,
                i + len(batch),
                exc,
            )
            # Continue with next batch rather than aborting entirely
            continue

    return results


def _property_to_dict(prop: dict[str, Any]) -> dict:
    """Convert a PubChem property record to the pipeline standard dict."""
    cid = prop.get("CID", "")
    return {
        "id": f"CID{cid}" if cid else "",
        "smiles": prop.get("CanonicalSMILES", ""),
        "name": prop.get("IUPACName", ""),
        "source": "PubChem",
        "molecular_weight": prop.get("MolecularWeight"),
        "activity_type": "active",
        "activity_value": None,
    }
