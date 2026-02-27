"""DrugBank open data parser.

DrugBank's full API requires a commercial licence and registered credentials.
This module supports two fallback strategies in order of preference:

  1. Local CSV file — downloaded from the DrugBank open data release at
     https://go.drugbank.com/releases/latest (free registration required).
     Expected path: ``{data_dir}/drugbank_open_structures.csv``
     Columns expected: drugbank_id, name, smiles, description

  2. DrugBank website scraping — a best-effort HTML scrape of the public
     search endpoint. Fragile; may break on layout changes.

Access requirements
-------------------
Full bioactivity and MoA data requires the DrugBank Academic or Commercial
licence. The open CSV only contains approved drug structures and basic metadata.
"""

from __future__ import annotations

import csv
import io
import logging
import os
from pathlib import Path

from drugdiscovery.utils.web import get_with_retries

logger = logging.getLogger(__name__)

_DRUGBANK_SEARCH_URL = "https://www.drugbank.ca/unearth/q?searcher=drugs&query={query}"


def search_drugbank(
    query: str,
    data_dir: str = "",
    limit: int = 100,
) -> list[dict]:
    """Search DrugBank open data for drugs matching a query string.

    Tries local CSV first, then falls back to a lightweight website search.
    Returns an empty list with a warning if neither source is available.

    Args:
        query: Target name or drug name to search for.
        data_dir: Directory containing ``drugbank_open_structures.csv``.
                  Defaults to the ``DRUGBANK_DATA_DIR`` environment variable,
                  then the current working directory.
        limit: Maximum number of results.

    Returns:
        List of dicts with keys:
            id, name, smiles, description, moa, source
    """
    data_dir = data_dir or os.environ.get("DRUGBANK_DATA_DIR", "")

    # Strategy 1: local CSV
    results = _search_local_csv(query, data_dir, limit)
    if results:
        logger.info("DrugBank search_drugbank(%r) from CSV: %d results", query, len(results))
        return results

    # Strategy 2: public website
    results = _search_website(query, limit)
    if results is not None:
        logger.info(
            "DrugBank search_drugbank(%r) from website: %d results", query, len(results)
        )
        return results

    logger.warning(
        "DrugBank: no local CSV found at %r and website search unavailable. "
        "Download the open data CSV from https://go.drugbank.com/releases/latest "
        "and set DRUGBANK_DATA_DIR or pass data_dir=...",
        data_dir or "<cwd>",
    )
    return []


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _search_local_csv(query: str, data_dir: str, limit: int) -> list[dict]:
    """Search the DrugBank open CSV for rows matching query (case-insensitive)."""
    csv_path = Path(data_dir) / "drugbank_open_structures.csv" if data_dir else Path("drugbank_open_structures.csv")
    if not csv_path.exists():
        return []

    results: list[dict] = []
    query_lower = query.lower()

    try:
        with csv_path.open(newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                searchable = " ".join(row.values()).lower()
                if query_lower in searchable:
                    results.append(_csv_row_to_dict(row))
                    if len(results) >= limit:
                        break
    except Exception as exc:  # noqa: BLE001
        logger.warning("DrugBank: failed to read local CSV %s: %s", csv_path, exc)
        return []

    return results


def _search_website(query: str, limit: int) -> list[dict] | None:
    """Attempt a lightweight HTML scrape of the DrugBank public search page.

    Returns None if the request fails (caller will warn and return []).
    Returns [] if request succeeded but no results were parseable.

    Note: This scrape is intentionally minimal to avoid overloading the server
    and to limit fragility. Only the drugbank_id, name, and a partial
    description are extracted from the search result page.
    """
    try:
        url = _DRUGBANK_SEARCH_URL.format(query=query)
        resp = get_with_retries(url, attempts=2, headers={"Accept": "text/html"})
        html = resp.text
    except Exception as exc:  # noqa: BLE001
        logger.warning("DrugBank website search failed: %s", exc)
        return None

    # Minimal extraction: look for DrugBank IDs in href links (DB\d{5})
    import re

    ids_seen: set[str] = set()
    results: list[dict] = []
    for match in re.finditer(r'href="/drugs/(DB\d{5})"[^>]*>([^<]+)</a>', html):
        db_id = match.group(1)
        name = match.group(2).strip()
        if db_id in ids_seen:
            continue
        ids_seen.add(db_id)
        results.append(
            {
                "id": db_id,
                "name": name,
                "smiles": "",          # not available without API access
                "description": "",
                "moa": "",
                "source": "DrugBank",
            }
        )
        if len(results) >= limit:
            break

    return results


def _csv_row_to_dict(row: dict) -> dict:
    """Normalise a CSV row to the pipeline standard dict."""
    return {
        "id": row.get("drugbank_id") or row.get("DrugBank ID", ""),
        "name": row.get("name") or row.get("Name", ""),
        "smiles": row.get("smiles") or row.get("SMILES", ""),
        "description": row.get("description") or row.get("Description", ""),
        "moa": row.get("mechanism_of_action") or row.get("Mechanism of Action", ""),
        "source": "DrugBank",
    }
