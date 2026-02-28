"""BindingDB API client for quantitative binding affinity data.

BindingDB contains ~3.2 million experimentally determined binding data points
including Ki, Kd, IC50, and EC50 measurements for drug-target interactions.

API endpoints used:
  - getLigandsByUniprots: search by UniProt accession
  - getLigandsByTarget: search by target name

No authentication is required for the public REST API.

Reference: https://www.bindingdb.org/rwd/bind/ByUniprotAccession.jsp
"""

from __future__ import annotations

import logging
from typing import Any

from drugdiscovery.utils.web import get_with_retries

logger = logging.getLogger(__name__)

BINDINGDB_BASE = "https://bindingdb.org/axis2/services/BDBService"

_TIMEOUT = 60  # BindingDB can be slow


def search_by_uniprot(uniprot_id: str, limit: int = 200) -> list[dict]:
    """Search BindingDB for ligands with measured binding to a UniProt target.

    Args:
        uniprot_id: UniProt accession, e.g. "P00533".
        limit: Maximum number of results to return.

    Returns:
        List of dicts with keys:
            id, smiles, name, source, ki, kd, ic50, ec50,
            affinity_type, affinity_value
        Returns an empty list on any error.
    """
    if not uniprot_id:
        return []

    try:
        url = f"{BINDINGDB_BASE}/getLigandsByUniprots"
        params: dict[str, Any] = {
            "uniprot": uniprot_id,
            "response": "json",
        }
        resp = get_with_retries(url, params=params, timeout=_TIMEOUT)
        data = resp.json()

        results = _parse_binding_response(data, limit)
        logger.info(
            "BindingDB search_by_uniprot(%r): %d results", uniprot_id, len(results)
        )
        return results

    except Exception as exc:  # noqa: BLE001
        logger.warning("BindingDB search_by_uniprot failed: %s", exc)
        return []


def search_by_target(target_name: str, limit: int = 200) -> list[dict]:
    """Search BindingDB for ligands with measured binding to a named target.

    Args:
        target_name: Human-readable target name, e.g. "EGFR".
        limit: Maximum number of results to return.

    Returns:
        List of dicts with keys:
            id, smiles, name, source, ki, kd, ic50, ec50,
            affinity_type, affinity_value
        Returns an empty list on any error.
    """
    if not target_name:
        return []

    try:
        url = f"{BINDINGDB_BASE}/getLigandsByTarget"
        params: dict[str, Any] = {
            "target": target_name,
            "response": "json",
        }
        resp = get_with_retries(url, params=params, timeout=_TIMEOUT)
        data = resp.json()

        results = _parse_binding_response(data, limit)
        logger.info(
            "BindingDB search_by_target(%r): %d results", target_name, len(results)
        )
        return results

    except Exception as exc:  # noqa: BLE001
        logger.warning("BindingDB search_by_target failed: %s", exc)
        return []


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _parse_binding_response(data: dict | list, limit: int) -> list[dict]:
    """Normalise a BindingDB API response to the pipeline standard list.

    The BindingDB JSON response wraps records in various structures depending
    on the endpoint and number of results.  This function handles:
      - ``{"getLigandsByUniprotsResponse": {"affinities": [...]}}``
      - ``{"getLigandsByTargetResponse": {"affinities": [...]}}``
      - A bare list of affinity records
      - A single record returned as a dict (instead of a one-element list)
    """
    # Unwrap known response envelopes
    if isinstance(data, dict):
        for key in (
            "getLigandsByUniprotsResponse",
            "getLigandsByTargetResponse",
        ):
            if key in data:
                data = data[key]
                break

        # At this level the records should be under "affinities"
        if isinstance(data, dict):
            records = data.get("affinities", data.get("data", []))
        else:
            records = data
    else:
        records = data

    # A single record may be returned as a dict instead of a list
    if isinstance(records, dict):
        records = [records]

    results: list[dict] = []
    for rec in records[:limit]:
        results.append(_record_to_dict(rec))
    return results


def _record_to_dict(rec: dict) -> dict:
    """Convert a raw BindingDB affinity record to the pipeline standard dict."""
    ki = _safe_affinity(rec.get("ki"))
    kd = _safe_affinity(rec.get("kd"))
    ic50 = _safe_affinity(rec.get("ic50"))
    ec50 = _safe_affinity(rec.get("ec50"))

    # Determine the best (lowest nM) affinity among available measures
    affinity_type, affinity_value = _best_affinity(ki, kd, ic50, ec50)

    identifier = (
        rec.get("monomerid")
        or rec.get("ligand_id")
        or rec.get("zinc_id")
        or rec.get("cid")
        or ""
    )

    return {
        "id": str(identifier),
        "smiles": rec.get("smiles") or rec.get("ligand_smiles", ""),
        "name": rec.get("ligand_name") or rec.get("name", ""),
        "source": "BindingDB",
        "ki": ki,
        "kd": kd,
        "ic50": ic50,
        "ec50": ec50,
        "affinity_type": affinity_type,
        "affinity_value": affinity_value,
    }


def _safe_affinity(val: Any) -> float | None:
    """Safely parse an affinity value to a float in nM.

    BindingDB may return values as strings, numbers, or entries like
    ">10000", "<1", or empty strings.  Relational prefixes are stripped
    and the numeric portion is returned.

    Returns:
        Parsed float value in nM, or None if parsing fails.
    """
    if val is None:
        return None

    val_str = str(val).strip()
    if not val_str:
        return None

    # Strip relational prefixes (e.g. ">", "<", ">=", "<=", "~")
    for prefix in (">=", "<=", ">", "<", "~"):
        if val_str.startswith(prefix):
            val_str = val_str[len(prefix):].strip()
            break

    try:
        parsed = float(val_str)
        return parsed if parsed >= 0 else None
    except (ValueError, TypeError):
        return None


def _best_affinity(
    ki: float | None,
    kd: float | None,
    ic50: float | None,
    ec50: float | None,
) -> tuple[str, float | None]:
    """Pick the strongest (lowest nM) affinity among available measures.

    Returns:
        Tuple of (affinity_type, affinity_value).  If no valid affinity is
        available, returns ("", None).
    """
    candidates: list[tuple[str, float]] = []
    if ki is not None:
        candidates.append(("Ki", ki))
    if kd is not None:
        candidates.append(("Kd", kd))
    if ic50 is not None:
        candidates.append(("IC50", ic50))
    if ec50 is not None:
        candidates.append(("EC50", ec50))

    if not candidates:
        return ("", None)

    best = min(candidates, key=lambda x: x[1])
    return best
