"""PubChem IUPAC name resolution for small-molecule candidates."""

from __future__ import annotations

import logging
import time
from urllib.parse import quote

from drugdiscovery.types import Candidate
from drugdiscovery.utils.web import get_with_retries

logger = logging.getLogger(__name__)

PUBCHEM_IUPAC_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles"
    "/{smiles}/property/IUPACName/JSON"
)

# PubChem rate limit: 5 requests/second; 0.2s delay stays within budget.
_RATE_DELAY = 0.2


def get_iupac_name(smiles: str) -> str:
    """Look up the IUPAC name for a SMILES string via PubChem REST API.

    Returns the IUPAC name string, or ``""`` on any failure.
    """
    if not smiles:
        return ""
    url = PUBCHEM_IUPAC_URL.format(smiles=quote(smiles, safe=""))
    try:
        resp = get_with_retries(url, timeout=15, attempts=2, delay=2)
        data = resp.json()
        props = data.get("PropertyTable", {}).get("Properties", [])
        if props:
            return props[0].get("IUPACName", "")
    except Exception:
        logger.debug("PubChem IUPAC lookup failed for %.40s", smiles)
    return ""


def resolve_iupac_names(candidates: list[Candidate]) -> list[Candidate]:
    """Batch-resolve IUPAC names for all small-molecule candidates.

    Peptide candidates are skipped (they use sequence notation).
    Results are cached per SMILES to avoid duplicate API calls for
    modified candidates that share a parent SMILES.

    Mutates candidates in-place and also returns the list for convenience.
    """
    cache: dict[str, str] = {}
    resolved = 0

    for c in candidates:
        if c.modality == "peptide" or not c.smiles:
            continue
        if c.iupac_name:
            continue  # already set

        if c.smiles in cache:
            c.iupac_name = cache[c.smiles]
            continue

        name = get_iupac_name(c.smiles)
        cache[c.smiles] = name
        c.iupac_name = name
        if name:
            resolved += 1
        time.sleep(_RATE_DELAY)

    logger.info("Resolved IUPAC names for %d / %d small-molecule candidates",
                resolved, sum(1 for c in candidates if c.modality != "peptide"))
    return candidates
