"""DrugBank database client.

Supports multiple data source strategies in order of preference:

  1. Pre-built index — a JSON index built from the full XML database for
     fast repeated searches. Created automatically on first XML parse.
     Path: ``{data_dir}/drugbank_index.json``

  2. Full XML database — the comprehensive DrugBank XML file containing
     all drug data including SMILES, targets, and mechanism of action.
     Download via (requires free academic account):
       curl -Lfv -o drugbank_full.zip -u EMAIL:PASSWORD \\
         https://go.drugbank.com/releases/5-1-14/downloads/all-full-database
       unzip drugbank_full.zip -d {data_dir}/
     Expected path: ``{data_dir}/full database.xml`` or ``{data_dir}/drugbank.xml``

  3. Open structures CSV — the free DrugBank open data CSV.
     Expected path: ``{data_dir}/drugbank_open_structures.csv``

  4. DrugBank website scraping — best-effort HTML scrape (fragile).

Environment:
  DRUGBANK_DATA_DIR — directory containing DrugBank data files.
"""

from __future__ import annotations

import csv
import json
import logging
import os
import re
import zipfile
from pathlib import Path
from typing import Any
from xml.etree.ElementTree import iterparse

from drugdiscovery.utils.web import get_with_retries

logger = logging.getLogger(__name__)

_DRUGBANK_SEARCH_URL = "https://www.drugbank.ca/unearth/q?searcher=drugs&query={query}"

# XML namespace used in DrugBank full database
_NS = "{http://www.drugbank.ca}"

# Default location for local data
_DEFAULT_DATA_DIR = Path("data/drugbank")


def search_drugbank(
    query: str,
    data_dir: str = "",
    limit: int = 100,
) -> list[dict]:
    """Search DrugBank for drugs matching a query string (target or drug name).

    Tries data sources in order: index -> XML -> CSV -> website.

    Args:
        query: Target gene name, drug name, or search term.
        data_dir: Directory containing DrugBank data files.
        limit: Maximum number of results.

    Returns:
        List of dicts with keys: id, name, smiles, description, moa, source
    """
    data_dir = _resolve_data_dir(data_dir)

    # Strategy 1: pre-built JSON index (fastest)
    results = _search_index(query, data_dir, limit)
    if results:
        logger.info("DrugBank search(%r) from index: %d results", query, len(results))
        return results

    # Strategy 2: full XML database (comprehensive, builds index for next time)
    results = _search_xml(query, data_dir, limit)
    if results:
        logger.info("DrugBank search(%r) from XML: %d results", query, len(results))
        return results

    # Strategy 3: open structures CSV
    results = _search_local_csv(query, data_dir, limit)
    if results:
        logger.info("DrugBank search(%r) from CSV: %d results", query, len(results))
        return results

    # Strategy 4: website scraping (fragile fallback)
    results = _search_website(query, limit)
    if results is not None:
        logger.info("DrugBank search(%r) from website: %d results", query, len(results))
        return results

    logger.warning(
        "DrugBank: no data sources found. Download the full database:\n"
        "  curl -Lfv -o drugbank_full.zip -u EMAIL:PASSWORD \\\n"
        "    https://go.drugbank.com/releases/5-1-14/downloads/all-full-database\n"
        "  unzip drugbank_full.zip -d %s/",
        data_dir or "data/drugbank",
    )
    return []


def build_index(data_dir: str = "") -> bool:
    """Build a JSON index from the full DrugBank XML for fast searching.

    Call this once after downloading and extracting the full database.

    Returns:
        True if index was built successfully.
    """
    data_dir = _resolve_data_dir(data_dir)
    xml_path = _find_xml(data_dir)
    if not xml_path:
        logger.error("No DrugBank XML found in %s", data_dir)
        return False

    logger.info("Building DrugBank index from %s (this may take a minute)...", xml_path)
    drugs = _parse_full_xml(xml_path)
    index_path = Path(data_dir) / "drugbank_index.json"
    with open(index_path, "w", encoding="utf-8") as f:
        json.dump(drugs, f, separators=(",", ":"))
    logger.info("DrugBank index built: %d drugs saved to %s", len(drugs), index_path)
    return True


# ---------------------------------------------------------------------------
# Data directory resolution
# ---------------------------------------------------------------------------


def _resolve_data_dir(data_dir: str) -> str:
    """Resolve the DrugBank data directory."""
    if data_dir:
        return data_dir

    env_dir = os.environ.get("DRUGBANK_DATA_DIR", "")
    if env_dir and Path(env_dir).is_dir():
        return env_dir

    if _DEFAULT_DATA_DIR.is_dir():
        return str(_DEFAULT_DATA_DIR)

    # Check relative to project root
    project_root = Path(__file__).resolve().parent.parent.parent
    candidate = project_root / "data" / "drugbank"
    if candidate.is_dir():
        return str(candidate)

    return ""


# ---------------------------------------------------------------------------
# Strategy 1: JSON index
# ---------------------------------------------------------------------------


def _search_index(query: str, data_dir: str, limit: int) -> list[dict]:
    """Search the pre-built JSON index."""
    if not data_dir:
        return []

    index_path = Path(data_dir) / "drugbank_index.json"
    if not index_path.is_file():
        return []

    try:
        with open(index_path, encoding="utf-8") as f:
            drugs: list[dict] = json.load(f)
    except Exception as exc:  # noqa: BLE001
        logger.warning("DrugBank: failed to read index: %s", exc)
        return []

    return _filter_drugs(drugs, query, limit)


# ---------------------------------------------------------------------------
# Strategy 2: Full XML database
# ---------------------------------------------------------------------------


def _find_xml(data_dir: str) -> Path | None:
    """Find the DrugBank XML file, handling ZIP extraction if needed."""
    if not data_dir:
        return None

    base = Path(data_dir)

    # Check for existing XML files
    xml_names = [
        "full database.xml",
        "drugbank.xml",
        "drugbank_all_full_database.xml",
    ]
    for name in xml_names:
        p = base / name
        if p.is_file():
            return p

    # Check for any .xml file
    for p in base.glob("*.xml"):
        if p.stat().st_size > 1_000_000:  # > 1MB likely full DB
            return p

    # Try to extract from ZIP
    zip_names = [
        "drugbank_full.zip",
        "drugbank.zip",
        "all-full-database.zip",
    ]
    for name in zip_names:
        zip_path = base / name
        if zip_path.is_file():
            return _extract_xml_from_zip(zip_path, base)

    # Check for any ZIP
    for p in base.glob("*.zip"):
        result = _extract_xml_from_zip(p, base)
        if result:
            return result

    return None


def _extract_xml_from_zip(zip_path: Path, dest_dir: Path) -> Path | None:
    """Extract the XML file from a DrugBank ZIP archive."""
    try:
        with zipfile.ZipFile(zip_path, "r") as zf:
            xml_files = [n for n in zf.namelist() if n.endswith(".xml")]
            if not xml_files:
                return None

            # Extract the largest XML (likely the full database)
            target = max(xml_files, key=lambda n: zf.getinfo(n).file_size)
            logger.info("Extracting DrugBank XML: %s", target)
            zf.extract(target, dest_dir)
            return dest_dir / target
    except Exception as exc:  # noqa: BLE001
        logger.warning("DrugBank: failed to extract ZIP %s: %s", zip_path, exc)
        return None


def _search_xml(query: str, data_dir: str, limit: int) -> list[dict]:
    """Search the full DrugBank XML and build index for future use."""
    xml_path = _find_xml(data_dir)
    if not xml_path:
        return []

    try:
        drugs = _parse_full_xml(xml_path)
    except Exception as exc:  # noqa: BLE001
        logger.warning("DrugBank: failed to parse XML %s: %s", xml_path, exc)
        return []

    # Build index for next time
    if data_dir:
        try:
            index_path = Path(data_dir) / "drugbank_index.json"
            with open(index_path, "w", encoding="utf-8") as f:
                json.dump(drugs, f, separators=(",", ":"))
            logger.info("DrugBank index saved: %d drugs", len(drugs))
        except Exception:  # noqa: BLE001
            pass

    return _filter_drugs(drugs, query, limit)


def _parse_full_xml(xml_path: Path) -> list[dict]:
    """Stream-parse the full DrugBank XML extracting drug records.

    Uses iterparse for memory-efficient processing of the large XML file.
    """
    drugs: list[dict] = []

    for event, elem in iterparse(str(xml_path), events=("end",)):
        if elem.tag != f"{_NS}drug" and elem.tag != "drug":
            continue

        # Only process top-level drug elements (not nested)
        drug = _extract_drug(elem)
        if drug:
            drugs.append(drug)

        # Free memory
        elem.clear()

    logger.info("DrugBank XML parsed: %d drugs extracted", len(drugs))
    return drugs


def _extract_drug(elem: Any) -> dict | None:
    """Extract a drug record from an XML element."""
    ns = _NS

    # DrugBank ID
    db_id = ""
    for id_elem in elem.findall(f"{ns}drugbank-id"):
        if id_elem.get("primary") == "true" or not db_id:
            db_id = (id_elem.text or "").strip()

    if not db_id:
        # Try without namespace
        for id_elem in elem.findall("drugbank-id"):
            if id_elem.get("primary") == "true" or not db_id:
                db_id = (id_elem.text or "").strip()

    if not db_id:
        return None

    name = _get_text(elem, f"{ns}name") or _get_text(elem, "name") or ""
    desc = _get_text(elem, f"{ns}description") or _get_text(elem, "description") or ""
    moa = (
        _get_text(elem, f"{ns}mechanism-of-action")
        or _get_text(elem, "mechanism-of-action")
        or ""
    )

    # SMILES from calculated-properties
    smiles = ""
    for prop_container in [f"{ns}calculated-properties", "calculated-properties"]:
        container = elem.find(prop_container)
        if container is not None:
            for prop in container:
                kind = _get_text(prop, f"{ns}kind") or _get_text(prop, "kind") or ""
                if kind == "SMILES":
                    smiles = (
                        _get_text(prop, f"{ns}value")
                        or _get_text(prop, "value")
                        or ""
                    )
                    break
        if smiles:
            break

    # Target gene names
    targets: list[str] = []
    for targets_container in [f"{ns}targets", "targets"]:
        container = elem.find(targets_container)
        if container is not None:
            for target in container:
                for poly_tag in [f"{ns}polypeptide", "polypeptide"]:
                    for poly in target.findall(poly_tag):
                        gene = (
                            _get_text(poly, f"{ns}gene-name")
                            or _get_text(poly, "gene-name")
                            or ""
                        )
                        if gene:
                            targets.append(gene.upper())
            break

    return {
        "id": db_id,
        "name": name,
        "smiles": smiles,
        "description": desc[:500] if desc else "",
        "moa": moa[:500] if moa else "",
        "source": "DrugBank",
        "targets": targets,
    }


def _get_text(elem: Any, tag: str) -> str | None:
    """Get text content of a child element."""
    child = elem.find(tag)
    if child is not None and child.text:
        return child.text.strip()
    return None


def _filter_drugs(drugs: list[dict], query: str, limit: int) -> list[dict]:
    """Filter drug list by query string (matches target genes, name, description)."""
    query_upper = query.upper().strip()
    results: list[dict] = []

    # First pass: exact target gene match (highest priority)
    for drug in drugs:
        if query_upper in drug.get("targets", []):
            results.append(_to_output(drug))
            if len(results) >= limit:
                return results

    # Second pass: name / description match
    for drug in drugs:
        if query_upper in drug.get("targets", []):
            continue  # already added
        searchable = f"{drug.get('name', '')} {drug.get('description', '')} {drug.get('moa', '')}".upper()
        if query_upper in searchable:
            results.append(_to_output(drug))
            if len(results) >= limit:
                return results

    return results


def _to_output(drug: dict) -> dict:
    """Convert internal drug dict to pipeline output format."""
    return {
        "id": drug.get("id", ""),
        "name": drug.get("name", ""),
        "smiles": drug.get("smiles", ""),
        "description": drug.get("description", ""),
        "moa": drug.get("moa", ""),
        "source": "DrugBank",
    }


# ---------------------------------------------------------------------------
# Strategy 3: Open structures CSV
# ---------------------------------------------------------------------------


def _search_local_csv(query: str, data_dir: str, limit: int) -> list[dict]:
    """Search the DrugBank open CSV for rows matching query."""
    candidates: list[Path] = []
    if data_dir:
        base = Path(data_dir)
        candidates.append(base / "drugbank_open_structures.csv")
        candidates.append(base / "drugbank_vocabulary.csv")
    candidates.append(Path("data") / "drugbank_open_structures.csv")
    candidates.append(Path("drugbank_open_structures.csv"))

    csv_path: Path | None = None
    for p in candidates:
        if p.exists():
            csv_path = p
            break

    if csv_path is None:
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
        logger.warning("DrugBank: failed to read CSV %s: %s", csv_path, exc)
        return []

    return results


# ---------------------------------------------------------------------------
# Strategy 4: Website scraping (fragile fallback)
# ---------------------------------------------------------------------------


def _search_website(query: str, limit: int) -> list[dict] | None:
    """Attempt a lightweight HTML scrape of the DrugBank search page."""
    try:
        url = _DRUGBANK_SEARCH_URL.format(query=query)
        resp = get_with_retries(url, attempts=2, headers={"Accept": "text/html"})
        html = resp.text
    except Exception as exc:  # noqa: BLE001
        logger.warning("DrugBank website search failed: %s", exc)
        return None

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
                "smiles": "",
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
