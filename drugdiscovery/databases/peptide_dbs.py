"""Peptide database clients: BIOPEP-UWM, APD3, and DBAASP.

These three databases catalogue bioactive peptides with varying coverage:

  - BIOPEP-UWM  https://biochemia.uwm.edu.pl/biopep-uwm/
    Focuses on food-derived bioactive peptides and enzymatic hydrolysates.
    No public REST API; this module performs lightweight HTML scraping.

  - APD3  https://aps.unmc.edu/
    Antimicrobial Peptide Database (third generation).
    Provides a query form interface; REST-style URLs are available for
    individual entries but bulk search requires form submission.

  - DBAASP  https://dbaasp.org/
    Database of Antimicrobial Activity and Structure of Peptides.
    Offers a JSON search API under /search/peptides.

Access notes
------------
BIOPEP-UWM and APD3 do not expose stable public REST APIs.  Scraping is
fragile; updates to website layout will silently return empty lists.  Check
logs for warnings if results are unexpectedly empty.

DBAASP exposes a reasonable JSON API and is the most reliable of the three.
"""

from __future__ import annotations

import logging
import re
from typing import Any

from drugdiscovery.utils.web import get_with_retries, post_with_retries

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Public endpoints
# ---------------------------------------------------------------------------

_BIOPEP_BASE = "https://biochemia.uwm.edu.pl/biopep-uwm"
_APD3_BASE = "https://aps.unmc.edu"
_DBAASP_BASE = "https://dbaasp.org"

_JSON_HEADERS = {"Accept": "application/json", "Content-Type": "application/json"}
_HTML_HEADERS = {"Accept": "text/html"}


# ===========================================================================
# BIOPEP-UWM
# ===========================================================================


def search_biopep(activity: str = "inhibitor", limit: int = 100) -> list[dict]:
    """Search BIOPEP-UWM for bioactive peptides by activity keyword.

    Because BIOPEP-UWM has no public REST API, this function scrapes the
    public search results page.  Results are limited in detail (no SMILES,
    since peptides are represented as sequences).

    Args:
        activity: Activity keyword, e.g. "inhibitor", "antihypertensive",
                  "antioxidant".  Sent as a URL query parameter.
        limit: Maximum number of results to return.

    Returns:
        List of dicts with keys:
            id, sequence, name, source, activity, length
        Returns an empty list on any error or if the site is unreachable.
    """
    try:
        url = f"{_BIOPEP_BASE}/search/"
        params: dict[str, Any] = {"activity": activity, "format": "json"}
        # Try JSON endpoint first (may exist in future versions)
        try:
            resp = get_with_retries(url, headers=_JSON_HEADERS, params=params, attempts=2)
            data = resp.json()
            results = _parse_biopep_json(data, activity, limit)
            if results:
                logger.info("BIOPEP search_biopep(%r): %d results", activity, len(results))
                return results
        except Exception:
            pass  # Fall through to HTML scrape

        # Fallback: HTML scrape
        results = _scrape_biopep_html(activity, limit)
        logger.info("BIOPEP search_biopep(%r) via HTML: %d results", activity, len(results))
        return results

    except Exception as exc:  # noqa: BLE001
        logger.warning("BIOPEP search_biopep failed: %s", exc)
        return []


# ===========================================================================
# APD3
# ===========================================================================


def search_apd3(activity: str = "antibacterial", limit: int = 100) -> list[dict]:
    """Search APD3 (Antimicrobial Peptide Database) for peptides by activity.

    APD3 exposes a PHP-based query form.  This function submits GET requests
    to the search endpoint and parses the HTML table response.

    Args:
        activity: Activity descriptor, e.g. "antibacterial", "antifungal",
                  "antiviral", "anticancer".
        limit: Maximum number of results.

    Returns:
        List of dicts with keys:
            id, sequence, name, source, activity, length
        Returns an empty list on any error.
    """
    try:
        url = f"{_APD3_BASE}/ap/info/dbsearch.php"
        params: dict[str, Any] = {
            "activity": activity,
            "b": "1",
            "end": str(limit),
        }
        resp = get_with_retries(url, headers=_HTML_HEADERS, params=params, attempts=2)
        results = _parse_apd3_html(resp.text, activity, limit)
        logger.info("APD3 search_apd3(%r): %d results", activity, len(results))
        return results

    except Exception as exc:  # noqa: BLE001
        logger.warning("APD3 search_apd3 failed: %s", exc)
        return []


# ===========================================================================
# DBAASP
# ===========================================================================


def search_dbaasp(
    activity: str = "",
    target: str = "",
    limit: int = 100,
) -> list[dict]:
    """Search DBAASP for peptides by activity and/or target organism.

    DBAASP exposes a JSON search API at /search/peptides (POST).
    At least one of ``activity`` or ``target`` should be provided.

    Args:
        activity: Activity filter, e.g. "antibacterial", "antifungal",
                  "anticancer".  Empty string omits this filter.
        target: Target organism or cell type, e.g. "Staphylococcus aureus",
                "cancer cell".  Empty string omits this filter.
        limit: Maximum number of results.

    Returns:
        List of dicts with keys:
            id, sequence, name, source, activity, length
        Returns an empty list on any error.
    """
    try:
        url = f"{_DBAASP_BASE}/search/peptides"
        payload: dict[str, Any] = {"page": 1, "pageSize": limit}
        if activity:
            payload["activity"] = activity
        if target:
            payload["target"] = target

        resp = post_with_retries(url, json=payload, headers=_JSON_HEADERS, attempts=2)
        data = resp.json()
        results = _parse_dbaasp_json(data, activity or target, limit)
        logger.info(
            "DBAASP search_dbaasp(activity=%r, target=%r): %d results",
            activity,
            target,
            len(results),
        )
        return results

    except Exception as exc:  # noqa: BLE001
        logger.warning("DBAASP search_dbaasp failed: %s", exc)
        return []


# ===========================================================================
# Internal parsers
# ===========================================================================


def _parse_biopep_json(data: dict | list, activity: str, limit: int) -> list[dict]:
    """Parse a BIOPEP JSON response (optimistic; format not publicly documented)."""
    if isinstance(data, list):
        records = data
    else:
        records = data.get("results") or data.get("peptides") or []

    results: list[dict] = []
    for rec in records[:limit]:
        seq = rec.get("sequence") or rec.get("seq", "")
        results.append(
            {
                "id": str(rec.get("id", "")),
                "sequence": seq,
                "name": rec.get("name") or rec.get("peptide_name", ""),
                "source": "BIOPEP-UWM",
                "activity": rec.get("activity") or activity,
                "length": len(seq),
            }
        )
    return results


def _scrape_biopep_html(activity: str, limit: int) -> list[dict]:
    """Scrape BIOPEP-UWM HTML search results for peptide sequences."""
    try:
        url = f"{_BIOPEP_BASE}/search/"
        params: dict[str, Any] = {"activity": activity}
        resp = get_with_retries(url, headers=_HTML_HEADERS, params=params, attempts=2)
        html = resp.text
    except Exception as exc:  # noqa: BLE001
        logger.debug("BIOPEP HTML scrape request failed: %s", exc)
        return []

    results: list[dict] = []
    # Heuristic: look for peptide sequence cells (uppercase letters, 2-50 chars)
    # and adjacent identifiers. Layout-dependent; may need updating.
    seq_pattern = re.compile(r"\b([ACDEFGHIKLMNPQRSTVWY]{2,50})\b")
    id_pattern = re.compile(r"BIOPEP[- ]?(\d+)", re.IGNORECASE)

    ids_found = id_pattern.findall(html)
    seqs_found = seq_pattern.findall(html)

    for i, seq in enumerate(seqs_found[:limit]):
        pep_id = f"BIOPEP-{ids_found[i]}" if i < len(ids_found) else f"BIOPEP-unknown-{i+1}"
        results.append(
            {
                "id": pep_id,
                "sequence": seq,
                "name": pep_id,
                "source": "BIOPEP-UWM",
                "activity": activity,
                "length": len(seq),
            }
        )

    return results


def _parse_apd3_html(html: str, activity: str, limit: int) -> list[dict]:
    """Extract peptide entries from APD3 HTML search result table."""
    results: list[dict] = []
    # APD3 table rows contain: AP ID | Name | Sequence | ... activity ...
    # Pattern: APD ID is "AP" followed by digits
    row_pattern = re.compile(
        r"AP(\d+)[^\w].*?"  # APD ID
        r"([ACDEFGHIKLMNPQRSTVWY]{4,})",  # sequence (amino acid letters)
        re.DOTALL,
    )
    for match in row_pattern.finditer(html):
        apd_id = f"APD{match.group(1)}"
        seq = match.group(2)
        results.append(
            {
                "id": apd_id,
                "sequence": seq,
                "name": apd_id,
                "source": "APD3",
                "activity": activity,
                "length": len(seq),
            }
        )
        if len(results) >= limit:
            break

    return results


def _parse_dbaasp_json(data: dict | list, activity_label: str, limit: int) -> list[dict]:
    """Normalise a DBAASP JSON response to the pipeline standard list."""
    if isinstance(data, list):
        records = data
    else:
        records = (
            data.get("peptides")
            or data.get("results")
            or data.get("data")
            or []
        )

    results: list[dict] = []
    for rec in records[:limit]:
        seq = rec.get("sequence") or rec.get("aa_sequence", "")
        # DBAASP returns a list of activities per peptide
        raw_activities = rec.get("activities") or rec.get("activity", [])
        if isinstance(raw_activities, list):
            activity_str = "; ".join(
                str(a.get("name", a) if isinstance(a, dict) else a)
                for a in raw_activities
            )
        else:
            activity_str = str(raw_activities) if raw_activities else activity_label

        results.append(
            {
                "id": str(rec.get("id") or rec.get("peptide_id", "")),
                "sequence": seq,
                "name": rec.get("name") or rec.get("peptide_name", ""),
                "source": "DBAASP",
                "activity": activity_str or activity_label,
                "length": len(seq),
            }
        )
    return results
