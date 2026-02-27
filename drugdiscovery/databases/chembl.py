"""ChEMBL REST API client for bioactivity data."""

from __future__ import annotations

import logging
from typing import Any

from drugdiscovery.utils.web import get_with_retries

logger = logging.getLogger(__name__)

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"

_HEADERS = {"Accept": "application/json"}


def search_by_target(
    target_name: str,
    organism: str = "Homo sapiens",
    limit: int = 500,
) -> list[dict]:
    """Search ChEMBL for compounds active against a named target.

    Workflow:
      1. Search for the target by name to obtain its ChEMBL target ID.
      2. Fetch bioactivity records filtered by pChEMBL >= 5 (micromolar or better).
      3. Collate molecule metadata and return a uniform list of dicts.

    Args:
        target_name: Human-readable target name, e.g. "YARS2".
        organism: Species filter applied to target search.
        limit: Maximum number of activity records to retrieve.

    Returns:
        List of dicts with keys:
            id, smiles, name, source, activity_type, activity_value, moa
        Returns an empty list on any error.
    """
    try:
        target_chembl_id = _resolve_target_id(target_name, organism)
        if not target_chembl_id:
            logger.warning(
                "ChEMBL: no target found for %r (organism=%s)", target_name, organism
            )
            return []

        activities = _fetch_activities(target_chembl_id, limit)
        results = _format_activities(activities)
        logger.info(
            "ChEMBL search_by_target(%r): %d results", target_name, len(results)
        )
        return results

    except Exception as exc:  # noqa: BLE001
        logger.warning("ChEMBL search_by_target failed: %s", exc)
        return []


def search_by_similarity(
    smiles: str,
    threshold: float = 70,
    limit: int = 100,
) -> list[dict]:
    """Search ChEMBL by Tanimoto similarity to a query SMILES.

    Uses the ChEMBL similarity endpoint which accepts a threshold in the range
    40-100 (percentage).

    Args:
        smiles: Query SMILES string.
        threshold: Minimum Tanimoto similarity (0-100). ChEMBL enforces >= 40.
        limit: Maximum number of results.

    Returns:
        List of dicts with keys:
            id, smiles, name, source, activity_type, activity_value, moa
        Returns an empty list on any error.
    """
    try:
        threshold_int = max(40, min(100, int(threshold)))
        url = f"{CHEMBL_BASE}/similarity/{smiles}/{threshold_int}.json"
        params: dict[str, Any] = {"limit": limit}

        resp = get_with_retries(url, headers=_HEADERS, params=params)
        data = resp.json()

        molecules = data.get("molecules", [])
        results = []
        for mol in molecules:
            results.append(_molecule_to_dict(mol))

        logger.info(
            "ChEMBL search_by_similarity(threshold=%s): %d results",
            threshold_int,
            len(results),
        )
        return results

    except Exception as exc:  # noqa: BLE001
        logger.warning("ChEMBL search_by_similarity failed: %s", exc)
        return []


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _resolve_target_id(target_name: str, organism: str) -> str | None:
    """Return the first ChEMBL target ID matching name + organism, or None."""
    url = f"{CHEMBL_BASE}/target/search.json"
    params: dict[str, Any] = {"q": target_name, "limit": 5}
    resp = get_with_retries(url, headers=_HEADERS, params=params)
    data = resp.json()

    targets = data.get("targets", [])
    for t in targets:
        if organism.lower() in t.get("organism", "").lower():
            return t.get("target_chembl_id")

    # Fallback: return the first result if no organism match
    if targets:
        logger.debug(
            "ChEMBL: organism filter %r did not match; using first target", organism
        )
        return targets[0].get("target_chembl_id")

    return None


def _fetch_activities(target_chembl_id: str, limit: int) -> list[dict]:
    """Fetch bioactivity records for a ChEMBL target ID."""
    url = f"{CHEMBL_BASE}/activity.json"
    params: dict[str, Any] = {
        "target_chembl_id": target_chembl_id,
        "pchembl_value__gte": 5,
        "limit": limit,
    }
    resp = get_with_retries(url, headers=_HEADERS, params=params)
    data = resp.json()
    return data.get("activities", [])


def _format_activities(activities: list[dict]) -> list[dict]:
    """Convert raw ChEMBL activity records to the pipeline standard dict."""
    results = []
    for act in activities:
        results.append(
            {
                "id": act.get("molecule_chembl_id", ""),
                "smiles": act.get("canonical_smiles", ""),
                "name": act.get("molecule_pref_name") or act.get("molecule_chembl_id", ""),
                "source": "ChEMBL",
                "activity_type": act.get("standard_type", ""),
                "activity_value": act.get("pchembl_value"),
                "moa": act.get("action_type") or "",
            }
        )
    return results


def _molecule_to_dict(mol: dict) -> dict:
    """Convert a ChEMBL molecule record (from similarity endpoint) to pipeline dict."""
    structures = mol.get("molecule_structures") or {}
    return {
        "id": mol.get("molecule_chembl_id", ""),
        "smiles": structures.get("canonical_smiles", ""),
        "name": mol.get("pref_name") or mol.get("molecule_chembl_id", ""),
        "source": "ChEMBL",
        "activity_type": "",
        "activity_value": None,
        "moa": "",
    }
