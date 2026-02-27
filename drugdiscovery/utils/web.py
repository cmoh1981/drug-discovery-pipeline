"""HTTP client with retry logic and common API helpers.

Refactored from YARS2 pipeline: fetch_sequence.py
"""

from __future__ import annotations

import logging
import time
from typing import Any, Optional

import requests
from requests.exceptions import ConnectionError, HTTPError, Timeout

logger = logging.getLogger(__name__)

DEFAULT_TIMEOUT = 60
DEFAULT_RETRIES = 3
DEFAULT_DELAY = 5


def get_with_retries(
    url: str,
    timeout: int = DEFAULT_TIMEOUT,
    attempts: int = DEFAULT_RETRIES,
    delay: int = DEFAULT_DELAY,
    headers: Optional[dict[str, str]] = None,
    params: Optional[dict[str, Any]] = None,
) -> requests.Response:
    """HTTP GET with retry logic and exponential backoff."""
    for attempt in range(1, attempts + 1):
        try:
            resp = requests.get(
                url, timeout=timeout, headers=headers, params=params
            )
            resp.raise_for_status()
            return resp
        except (HTTPError, ConnectionError, Timeout) as exc:
            logger.warning(
                "GET %s attempt %d/%d failed: %s", url, attempt, attempts, exc
            )
            if attempt < attempts:
                time.sleep(delay * attempt)  # exponential backoff
            else:
                raise
    raise RuntimeError(f"Failed to GET {url} after {attempts} attempts")


def post_with_retries(
    url: str,
    json: Optional[dict] = None,
    data: Optional[Any] = None,
    timeout: int = DEFAULT_TIMEOUT,
    attempts: int = DEFAULT_RETRIES,
    delay: int = DEFAULT_DELAY,
    headers: Optional[dict[str, str]] = None,
) -> requests.Response:
    """HTTP POST with retry logic."""
    for attempt in range(1, attempts + 1):
        try:
            resp = requests.post(
                url, json=json, data=data, timeout=timeout, headers=headers
            )
            resp.raise_for_status()
            return resp
        except (HTTPError, ConnectionError, Timeout) as exc:
            logger.warning(
                "POST %s attempt %d/%d failed: %s", url, attempt, attempts, exc
            )
            if attempt < attempts:
                time.sleep(delay * attempt)
            else:
                raise
    raise RuntimeError(f"Failed to POST {url} after {attempts} attempts")


# ---------------------------------------------------------------------------
# UniProt API helpers
# ---------------------------------------------------------------------------

UNIPROT_BASE = "https://rest.uniprot.org"
UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"
UNIPROT_JSON = "https://rest.uniprot.org/uniprotkb/{accession}.json"
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
ALPHAFOLD_PDB = "https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v4.pdb"


def fetch_uniprot_fasta(accession: str) -> str:
    """Fetch FASTA sequence from UniProt."""
    url = UNIPROT_FASTA.format(accession=accession)
    resp = get_with_retries(url)
    return resp.text


def fetch_uniprot_json(accession: str) -> dict:
    """Fetch full JSON record from UniProt."""
    url = UNIPROT_JSON.format(accession=accession)
    resp = get_with_retries(url)
    return resp.json()


def search_uniprot(query: str, limit: int = 5, organism: str = "9606") -> list[dict]:
    """Search UniProt and return matching entries.

    Args:
        query: Search term (gene name, protein name, etc.)
        limit: Max results to return
        organism: Taxonomy ID (9606 = human)
    """
    params = {
        "query": f"({query}) AND (organism_id:{organism})",
        "format": "json",
        "size": str(limit),
        "fields": "accession,gene_names,protein_name,organism_name,length",
    }
    resp = get_with_retries(UNIPROT_SEARCH, params=params)
    data = resp.json()
    return data.get("results", [])


def fetch_alphafold_pdb(accession: str) -> Optional[str]:
    """Fetch AlphaFold predicted structure PDB.

    Returns PDB text or None if not available.
    """
    url = ALPHAFOLD_PDB.format(accession=accession)
    try:
        resp = get_with_retries(url, attempts=2)
        return resp.text
    except Exception:
        logger.warning("AlphaFold structure not available for %s", accession)
        return None
