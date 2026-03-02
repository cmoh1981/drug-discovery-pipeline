"""BindingDB binding affinity data client.

BindingDB contains ~3.2 million experimentally determined binding data points
including Ki, Kd, IC50, and EC50 measurements for drug-target interactions.

Data sources (in priority order):
  1. Local TSV files from dhimmel/bindingdb (data/bindingdb/data/)
     - bindings-drugbank-collapsed.tsv: search by UniProt accession
     - bindings-drugbank-gene.tsv: search by gene symbol
  2. BindingDB REST API (fallback, may be unreliable)

SMILES are resolved via PubChem for local data hits (which only have
DrugBank IDs and drug names).

Setup:
  git clone --depth 1 https://github.com/dhimmel/bindingdb.git data/bindingdb
"""

from __future__ import annotations

import csv
import logging
import os
from pathlib import Path
from typing import Any

from drugdiscovery.utils.web import get_with_retries

logger = logging.getLogger(__name__)

BINDINGDB_BASE = "https://bindingdb.org/axis2/services/BDBService"

_TIMEOUT = 60

# Circuit breaker: after first API failure (e.g. expired SSL cert),
# skip all subsequent API calls in this process to avoid massive slowdowns.
_api_unavailable = False

# Default location for local data (relative to project root)
_DEFAULT_DATA_DIR = Path("data/bindingdb/data")


def search_by_uniprot(uniprot_id: str, limit: int = 200) -> list[dict]:
    """Search for ligands with measured binding to a UniProt target.

    Tries local TSV data first, then falls back to the BindingDB API.

    Args:
        uniprot_id: UniProt accession, e.g. "P00533".
        limit: Maximum number of results to return.

    Returns:
        List of dicts with keys:
            id, smiles, name, source, ki, kd, ic50, ec50,
            affinity_type, affinity_value
    """
    if not uniprot_id:
        return []

    # Strategy 1: local TSV data
    data_dir = _find_data_dir()
    if data_dir:
        results = _search_local_by_uniprot(data_dir, uniprot_id, limit)
        if results:
            # Resolve SMILES from PubChem for local hits
            _resolve_smiles_batch(results)
            logger.info(
                "BindingDB local search_by_uniprot(%r): %d results",
                uniprot_id,
                len(results),
            )
            return results

    # Strategy 2: BindingDB REST API (may be unreliable)
    return _api_search_by_uniprot(uniprot_id, limit)


def search_by_target(target_name: str, limit: int = 200) -> list[dict]:
    """Search for ligands with measured binding to a named target.

    Tries local TSV data first (gene symbol match), then API fallback.

    Args:
        target_name: Gene symbol or target name, e.g. "EGFR".
        limit: Maximum number of results to return.

    Returns:
        List of dicts with keys:
            id, smiles, name, source, ki, kd, ic50, ec50,
            affinity_type, affinity_value
    """
    if not target_name:
        return []

    # Strategy 1: local TSV data (gene symbol match)
    data_dir = _find_data_dir()
    if data_dir:
        results = _search_local_by_gene(data_dir, target_name, limit)
        if results:
            _resolve_smiles_batch(results)
            logger.info(
                "BindingDB local search_by_target(%r): %d results",
                target_name,
                len(results),
            )
            return results

    # Strategy 2: BindingDB REST API
    return _api_search_by_target(target_name, limit)


# ---------------------------------------------------------------------------
# Local TSV data search
# ---------------------------------------------------------------------------


def _find_data_dir() -> Path | None:
    """Locate the local bindingdb data directory."""
    # Check environment variable first
    env_dir = os.environ.get("BINDINGDB_DATA_DIR", "")
    if env_dir:
        p = Path(env_dir)
        if p.is_dir():
            return p

    # Check default location relative to CWD
    if _DEFAULT_DATA_DIR.is_dir():
        return _DEFAULT_DATA_DIR

    # Check relative to this file's location (project root)
    project_root = Path(__file__).resolve().parent.parent.parent
    candidate = project_root / "data" / "bindingdb" / "data"
    if candidate.is_dir():
        return candidate

    return None


def _search_local_by_uniprot(
    data_dir: Path, uniprot_id: str, limit: int
) -> list[dict]:
    """Search bindings-drugbank-collapsed.tsv by UniProt accession."""
    tsv_path = data_dir / "bindings-drugbank-collapsed.tsv"
    if not tsv_path.is_file():
        return []

    uniprot_upper = uniprot_id.upper().strip()
    results: list[dict] = []

    try:
        with tsv_path.open(encoding="utf-8", errors="replace") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                if row.get("uniprot", "").upper().strip() == uniprot_upper:
                    results.append(_collapsed_row_to_dict(row))
                    if len(results) >= limit:
                        break
    except Exception as exc:  # noqa: BLE001
        logger.warning("BindingDB: failed to read local TSV %s: %s", tsv_path, exc)
        return []

    return results


def _search_local_by_gene(
    data_dir: Path, gene_symbol: str, limit: int
) -> list[dict]:
    """Search bindings-drugbank-gene.tsv by gene symbol."""
    tsv_path = data_dir / "bindings-drugbank-gene.tsv"
    if not tsv_path.is_file():
        return []

    gene_upper = gene_symbol.upper().strip()
    results: list[dict] = []

    try:
        with tsv_path.open(encoding="utf-8", errors="replace") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                if row.get("gene_symbol", "").upper().strip() == gene_upper:
                    results.append(_gene_row_to_dict(row))
                    if len(results) >= limit:
                        break
    except Exception as exc:  # noqa: BLE001
        logger.warning("BindingDB: failed to read local TSV %s: %s", tsv_path, exc)
        return []

    return results


def _collapsed_row_to_dict(row: dict) -> dict:
    """Convert a collapsed TSV row to the pipeline standard dict."""
    measure = row.get("measure", "").strip()
    affinity = _safe_affinity(row.get("affinity_nM"))

    ki = affinity if measure == "Ki" else None
    kd = affinity if measure == "Kd" else None
    ic50 = affinity if measure == "IC50" else None
    ec50 = affinity if measure == "EC50" else None

    drugbank_id = row.get("drugbank_id", "")
    bindingdb_id = row.get("bindingdb_id", "")

    return {
        "id": drugbank_id or bindingdb_id,
        "smiles": "",  # resolved later via PubChem
        "name": drugbank_id,  # will be enriched if gene TSV has name
        "source": "BindingDB",
        "ki": ki,
        "kd": kd,
        "ic50": ic50,
        "ec50": ec50,
        "affinity_type": measure,
        "affinity_value": affinity,
        "drugbank_id": drugbank_id,
    }


def _gene_row_to_dict(row: dict) -> dict:
    """Convert a gene TSV row to the pipeline standard dict."""
    affinity = _safe_affinity(row.get("affinity_nM"))
    drugbank_id = row.get("drugbank_id", "")
    drug_name = row.get("drugbank_name", "")

    return {
        "id": drugbank_id,
        "smiles": "",  # resolved later via PubChem
        "name": drug_name or drugbank_id,
        "source": "BindingDB",
        "ki": affinity,  # gene TSV doesn't specify measure type
        "kd": None,
        "ic50": None,
        "ec50": None,
        "affinity_type": "affinity",
        "affinity_value": affinity,
        "drugbank_id": drugbank_id,
        "approved": row.get("drugbank_approved", "0") == "1",
    }


# ---------------------------------------------------------------------------
# SMILES resolution via PubChem
# ---------------------------------------------------------------------------


def _resolve_smiles_batch(results: list[dict]) -> None:
    """Resolve SMILES for results using PubChem name lookup.

    Modifies results in-place, setting the 'smiles' field.
    """
    names_to_resolve: list[tuple[int, str]] = []
    for i, r in enumerate(results):
        if not r.get("smiles"):
            name = r.get("name", "")
            if name and not name.startswith("DB"):
                names_to_resolve.append((i, name))

    if not names_to_resolve:
        # Try DrugBank IDs via PubChem
        for i, r in enumerate(results):
            if not r.get("smiles") and r.get("drugbank_id"):
                names_to_resolve.append((i, r["drugbank_id"]))

    resolved = 0
    for idx, name in names_to_resolve:
        smiles = _pubchem_smiles_lookup(name)
        if smiles:
            results[idx]["smiles"] = smiles
            resolved += 1

    if names_to_resolve:
        logger.info(
            "BindingDB: resolved SMILES for %d/%d compounds via PubChem",
            resolved,
            len(names_to_resolve),
        )


def _pubchem_smiles_lookup(name: str) -> str:
    """Look up canonical SMILES from PubChem by compound name or DrugBank ID."""
    if not name:
        return ""
    try:
        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
            f"{name}/property/CanonicalSMILES/JSON"
        )
        resp = get_with_retries(url, attempts=1, timeout=10, delay=1)
        data = resp.json()
        props = data.get("PropertyTable", {}).get("Properties", [])
        if props:
            p = props[0]
            return (
                p.get("CanonicalSMILES", "")
                or p.get("ConnectivitySMILES", "")
                or p.get("IsomericSMILES", "")
            )
    except Exception:  # noqa: BLE001
        pass
    return ""


# ---------------------------------------------------------------------------
# BindingDB REST API (fallback)
# ---------------------------------------------------------------------------


def _api_search_by_uniprot(uniprot_id: str, limit: int) -> list[dict]:
    """Search BindingDB REST API by UniProt accession (may be unreliable)."""
    global _api_unavailable  # noqa: PLW0603
    if _api_unavailable:
        return []

    try:
        url = f"{BINDINGDB_BASE}/getLigandsByUniprots"
        params: dict[str, Any] = {
            "uniprot": uniprot_id,
            "response": "json",
        }
        resp = get_with_retries(url, params=params, timeout=_TIMEOUT)
        data = resp.json()

        results = _parse_api_response(data, limit)
        logger.info(
            "BindingDB API search_by_uniprot(%r): %d results",
            uniprot_id,
            len(results),
        )
        return results

    except Exception as exc:  # noqa: BLE001
        _api_unavailable = True
        logger.warning(
            "BindingDB API unavailable (circuit breaker tripped, "
            "skipping future API calls): %s", exc,
        )
        return []


def _api_search_by_target(target_name: str, limit: int) -> list[dict]:
    """Search BindingDB REST API by target name (may be unreliable)."""
    global _api_unavailable  # noqa: PLW0603
    if _api_unavailable:
        return []

    try:
        url = f"{BINDINGDB_BASE}/getLigandsByTarget"
        params: dict[str, Any] = {
            "target": target_name,
            "response": "json",
        }
        resp = get_with_retries(url, params=params, timeout=_TIMEOUT)
        data = resp.json()

        results = _parse_api_response(data, limit)
        logger.info(
            "BindingDB API search_by_target(%r): %d results",
            target_name,
            len(results),
        )
        return results

    except Exception as exc:  # noqa: BLE001
        _api_unavailable = True
        logger.warning(
            "BindingDB API unavailable (circuit breaker tripped, "
            "skipping future API calls): %s", exc,
        )
        return []


def _parse_api_response(data: dict | list, limit: int) -> list[dict]:
    """Normalise a BindingDB API response to the pipeline standard list."""
    if isinstance(data, dict):
        for key in (
            "getLigandsByUniprotsResponse",
            "getLigandsByTargetResponse",
        ):
            if key in data:
                data = data[key]
                break

        if isinstance(data, dict):
            records = data.get("affinities", data.get("data", []))
        else:
            records = data
    else:
        records = data

    if isinstance(records, dict):
        records = [records]

    results: list[dict] = []
    for rec in records[:limit]:
        results.append(_api_record_to_dict(rec))
    return results


def _api_record_to_dict(rec: dict) -> dict:
    """Convert a raw BindingDB API affinity record to the pipeline standard dict."""
    ki = _safe_affinity(rec.get("ki"))
    kd = _safe_affinity(rec.get("kd"))
    ic50 = _safe_affinity(rec.get("ic50"))
    ec50 = _safe_affinity(rec.get("ec50"))

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


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _safe_affinity(val: Any) -> float | None:
    """Safely parse an affinity value to a float in nM."""
    if val is None:
        return None

    val_str = str(val).strip()
    if not val_str:
        return None

    for prefix in (">=", "<=", ">", "<", "~"):
        if val_str.startswith(prefix):
            val_str = val_str[len(prefix) :].strip()
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
    """Pick the strongest (lowest nM) affinity among available measures."""
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
