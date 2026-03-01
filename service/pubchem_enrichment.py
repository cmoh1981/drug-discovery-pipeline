"""PubChem compound enrichment for candidate details."""

from __future__ import annotations

import logging
import time
from urllib.parse import quote

import requests

logger = logging.getLogger(__name__)

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
_RATE_DELAY = 0.25
_TIMEOUT = 15

# Properties available from PubChem PUG REST
_PROPERTY_LIST = ",".join([
    "IUPACName",
    "MolecularFormula",
    "MolecularWeight",
    "CanonicalSMILES",
    "InChIKey",
    "XLogP",
    "ExactMass",
    "TPSA",
    "Complexity",
    "HBondDonorCount",
    "HBondAcceptorCount",
    "RotatableBondCount",
    "HeavyAtomCount",
])


def enrich_smiles(smiles: str) -> dict | None:
    """Look up a SMILES string in PubChem and return enriched data.

    Returns None if the compound is not found.
    """
    if not smiles or not smiles.strip():
        return None

    cid = _resolve_cid(smiles)
    if not cid:
        return None

    properties = _get_properties(cid)
    synonyms = _get_synonyms(cid)

    result = {
        "cid": cid,
        "pubchem_url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}",
        "iupac_name": "",
        "molecular_formula": "",
        "molecular_weight": None,
        "canonical_smiles": "",
        "inchikey": "",
        "xlogp": None,
        "exact_mass": None,
        "tpsa": None,
        "complexity": None,
        "h_bond_donor_count": None,
        "h_bond_acceptor_count": None,
        "rotatable_bond_count": None,
        "heavy_atom_count": None,
        "synonyms": [],
    }

    if properties:
        result["iupac_name"] = properties.get("IUPACName", "")
        result["molecular_formula"] = properties.get("MolecularFormula", "")
        result["canonical_smiles"] = properties.get("CanonicalSMILES", properties.get("ConnectivitySMILES", ""))
        result["inchikey"] = properties.get("InChIKey", "")
        result["xlogp"] = _safe_float(properties.get("XLogP"))
        result["tpsa"] = _safe_float(properties.get("TPSA"))
        result["complexity"] = _safe_float(properties.get("Complexity"))
        result["h_bond_donor_count"] = _safe_int(properties.get("HBondDonorCount"))
        result["h_bond_acceptor_count"] = _safe_int(properties.get("HBondAcceptorCount"))
        result["rotatable_bond_count"] = _safe_int(properties.get("RotatableBondCount"))
        result["heavy_atom_count"] = _safe_int(properties.get("HeavyAtomCount"))
        result["molecular_weight"] = _safe_float(properties.get("MolecularWeight"))
        result["exact_mass"] = _safe_float(properties.get("ExactMass"))

    if synonyms:
        result["synonyms"] = synonyms[:20]  # Cap at 20

    # Lipinski's Rule of Five assessment
    result["lipinski"] = _assess_lipinski(result)

    return result


def enrich_batch(smiles_list: list[str]) -> dict[str, dict | None]:
    """Enrich multiple SMILES strings. Returns {smiles: enrichment_data}."""
    results: dict[str, dict | None] = {}
    cache: dict[str, dict | None] = {}

    for smiles in smiles_list:
        if smiles in cache:
            results[smiles] = cache[smiles]
            continue
        data = enrich_smiles(smiles)
        cache[smiles] = data
        results[smiles] = data
        time.sleep(_RATE_DELAY)

    return results


def _resolve_cid(smiles: str) -> int | None:
    """Resolve a SMILES string to a PubChem CID."""
    encoded = quote(smiles, safe="")
    url = f"{PUBCHEM_BASE}/compound/smiles/{encoded}/cids/JSON"
    try:
        resp = requests.get(url, timeout=_TIMEOUT)
        if resp.status_code != 200:
            return None
        data = resp.json()
        cids = data.get("IdentifierList", {}).get("CID", [])
        return cids[0] if cids else None
    except Exception:
        logger.debug("PubChem CID lookup failed for %.40s", smiles)
        return None


def _get_properties(cid: int) -> dict | None:
    """Fetch compound properties from PubChem."""
    url = f"{PUBCHEM_BASE}/compound/cid/{cid}/property/{_PROPERTY_LIST}/JSON"
    time.sleep(_RATE_DELAY)
    try:
        resp = requests.get(url, timeout=_TIMEOUT)
        if resp.status_code != 200:
            return None
        data = resp.json()
        props = data.get("PropertyTable", {}).get("Properties", [])
        return props[0] if props else None
    except Exception:
        logger.debug("PubChem property fetch failed for CID %d", cid)
        return None


def _get_synonyms(cid: int) -> list[str]:
    """Fetch compound synonyms (common names, drug names, registry numbers)."""
    url = f"{PUBCHEM_BASE}/compound/cid/{cid}/synonyms/JSON"
    time.sleep(_RATE_DELAY)
    try:
        resp = requests.get(url, timeout=_TIMEOUT)
        if resp.status_code != 200:
            return []
        data = resp.json()
        info = data.get("InformationList", {}).get("Information", [])
        return info[0].get("Synonym", []) if info else []
    except Exception:
        logger.debug("PubChem synonym fetch failed for CID %d", cid)
        return []


def _assess_lipinski(data: dict) -> dict:
    """Assess Lipinski's Rule of Five for oral drug-likeness."""
    mw = data.get("molecular_weight")
    xlogp = data.get("xlogp")
    hbd = data.get("h_bond_donor_count")
    hba = data.get("h_bond_acceptor_count")

    violations = 0
    rules = {}

    if mw is not None:
        passed = mw <= 500
        rules["mw_le_500"] = {"value": mw, "passed": passed}
        if not passed:
            violations += 1

    if xlogp is not None:
        passed = xlogp <= 5
        rules["xlogp_le_5"] = {"value": xlogp, "passed": passed}
        if not passed:
            violations += 1

    if hbd is not None:
        passed = hbd <= 5
        rules["hbd_le_5"] = {"value": hbd, "passed": passed}
        if not passed:
            violations += 1

    if hba is not None:
        passed = hba <= 10
        rules["hba_le_10"] = {"value": hba, "passed": passed}
        if not passed:
            violations += 1

    return {
        "violations": violations,
        "drug_like": violations <= 1,
        "rules": rules,
    }


def _safe_float(value) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (ValueError, TypeError):
        return None


def _safe_int(value) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except (ValueError, TypeError):
        return None
