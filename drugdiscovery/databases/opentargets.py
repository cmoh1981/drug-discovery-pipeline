"""Open Targets Platform GraphQL API client.

Provides target-disease associations, known drugs, and clinical evidence
from the Open Targets Platform (https://platform.opentargets.org/).

All queries use the public GraphQL endpoint (no authentication required).
"""

from __future__ import annotations

import logging
from typing import Any

from drugdiscovery.utils.web import post_with_retries

logger = logging.getLogger(__name__)

OPENTARGETS_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"

_HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}

# ---------------------------------------------------------------------------
# GraphQL query templates
# ---------------------------------------------------------------------------

_TARGET_SEARCH_QUERY = """\
query TargetSearch($queryString: String!) {
  search(queryString: $queryString, entityNames: ["target"], page: {size: 5, index: 0}) {
    hits {
      id
    }
  }
}
"""

_ASSOCIATED_DISEASES_QUERY = """\
query AssociatedDiseases($ensemblId: String!, $size: Int!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    associatedDiseases(page: {size: $size, index: 0}) {
      rows {
        disease {
          id
          name
          therapeuticAreas {
            id
            name
          }
        }
        score
      }
    }
  }
}
"""

_KNOWN_DRUGS_QUERY = """\
query KnownDrugs($ensemblId: String!, $size: Int!) {
  target(ensemblId: $ensemblId) {
    knownDrugs(size: $size) {
      rows {
        drugId
        prefName
        drugType
        mechanismOfAction
        phase
      }
    }
  }
}
"""

_DISEASE_TARGETS_QUERY = """\
query DiseaseTargets($diseaseId: String!, $size: Int!) {
  disease(efoId: $diseaseId) {
    associatedTargets(page: {size: $size, index: 0}) {
      rows {
        target {
          approvedSymbol
        }
        score
      }
    }
  }
}
"""


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def get_target_diseases(gene_name: str, limit: int = 20) -> list[dict]:
    """Get diseases associated with a gene target from Open Targets.

    Args:
        gene_name: Human gene symbol, e.g. "EGFR".
        limit: Maximum number of disease associations to retrieve.

    Returns:
        List of dicts with keys:
            disease_id, disease_name, score, therapeutic_areas
        Returns an empty list on any error.
    """
    try:
        ensembl_id = _resolve_target_id(gene_name)
        if not ensembl_id:
            logger.warning(
                "OpenTargets: no target found for %r", gene_name
            )
            return []

        data = _graphql_request(
            _ASSOCIATED_DISEASES_QUERY,
            {"ensemblId": ensembl_id, "size": limit},
        )

        target = data.get("data", {}).get("target")
        if not target:
            logger.warning("OpenTargets: target query returned no data for %r", gene_name)
            return []

        rows = (
            target
            .get("associatedDiseases", {})
            .get("rows", [])
        )

        results = []
        for row in rows:
            disease = row.get("disease", {})
            therapeutic_areas = [
                area.get("name", "")
                for area in disease.get("therapeuticAreas", [])
            ]
            results.append(
                {
                    "disease_id": disease.get("id", ""),
                    "disease_name": disease.get("name", ""),
                    "score": row.get("score", 0.0),
                    "therapeutic_areas": therapeutic_areas,
                }
            )

        logger.info(
            "OpenTargets get_target_diseases(%r): %d results",
            gene_name,
            len(results),
        )
        return results

    except Exception as exc:  # noqa: BLE001
        logger.warning("OpenTargets get_target_diseases failed: %s", exc)
        return []


def get_known_drugs(gene_name: str, limit: int = 50) -> list[dict]:
    """Get known drugs and clinical compounds for a gene target.

    Args:
        gene_name: Human gene symbol, e.g. "EGFR".
        limit: Maximum number of drug records to retrieve.

    Returns:
        List of dicts with keys:
            id, drug_name, drug_type, phase, mechanism_of_action, source
        Returns an empty list on any error.
    """
    try:
        ensembl_id = _resolve_target_id(gene_name)
        if not ensembl_id:
            logger.warning(
                "OpenTargets: no target found for %r", gene_name
            )
            return []

        data = _graphql_request(
            _KNOWN_DRUGS_QUERY,
            {"ensemblId": ensembl_id, "size": limit},
        )

        target = data.get("data", {}).get("target")
        if not target:
            logger.warning("OpenTargets: knownDrugs query returned no data for %r", gene_name)
            return []

        rows = (
            target
            .get("knownDrugs", {})
            .get("rows", [])
        )

        results = []
        seen_ids: set[str] = set()
        for row in rows:
            drug_id = row.get("drugId", "")
            if drug_id in seen_ids:
                continue
            seen_ids.add(drug_id)

            results.append(
                {
                    "id": drug_id,
                    "drug_name": row.get("prefName") or drug_id,
                    "drug_type": row.get("drugType", ""),
                    "phase": row.get("phase"),
                    "mechanism_of_action": row.get("mechanismOfAction", ""),
                    "source": "OpenTargets",
                }
            )

        logger.info(
            "OpenTargets get_known_drugs(%r): %d unique drugs",
            gene_name,
            len(results),
        )
        return results

    except Exception as exc:  # noqa: BLE001
        logger.warning("OpenTargets get_known_drugs failed: %s", exc)
        return []


def get_disease_genes(gene_name: str, limit: int = 100) -> list[str]:
    """Get genes associated with diseases of a given target.

    Workflow:
      1. Resolve gene name to Ensembl ID.
      2. Fetch top diseases associated with the target.
      3. For each disease, fetch the associated gene targets.
      4. Collect and deduplicate gene symbols.

    Args:
        gene_name: Human gene symbol, e.g. "EGFR".
        limit: Maximum number of genes to collect per disease query.

    Returns:
        Deduplicated list of gene symbols. Returns an empty list on any error.
    """
    try:
        diseases = get_target_diseases(gene_name, limit=10)
        if not diseases:
            logger.info(
                "OpenTargets get_disease_genes(%r): no diseases found, returning empty",
                gene_name,
            )
            return []

        all_genes: set[str] = set()

        for disease in diseases[:5]:  # top 5 diseases to avoid excessive queries
            disease_id = disease.get("disease_id", "")
            if not disease_id:
                continue

            data = _graphql_request(
                _DISEASE_TARGETS_QUERY,
                {"diseaseId": disease_id, "size": limit},
            )

            disease_node = data.get("data", {}).get("disease")
            if not disease_node:
                continue

            rows = (
                disease_node
                .get("associatedTargets", {})
                .get("rows", [])
            )

            for row in rows:
                symbol = (
                    row.get("target", {})
                    .get("approvedSymbol", "")
                )
                if symbol:
                    all_genes.add(symbol)

        # Remove the query gene itself from the results
        all_genes.discard(gene_name)

        result = sorted(all_genes)
        logger.info(
            "OpenTargets get_disease_genes(%r): %d unique genes from %d diseases",
            gene_name,
            len(result),
            len(diseases),
        )
        return result

    except Exception as exc:  # noqa: BLE001
        logger.warning("OpenTargets get_disease_genes failed: %s", exc)
        return []


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _resolve_target_id(gene_name: str) -> str | None:
    """Map a gene symbol to an Ensembl gene ID (ENSG...) via Open Targets search.

    Args:
        gene_name: Human gene symbol, e.g. "EGFR".

    Returns:
        Ensembl gene ID string, or None if not found.
    """
    data = _graphql_request(
        _TARGET_SEARCH_QUERY,
        {"queryString": gene_name},
    )

    hits = (
        data.get("data", {})
        .get("search", {})
        .get("hits", [])
    )

    if not hits:
        return None

    target_id = hits[0].get("id", "")
    if target_id.startswith("ENSG"):
        logger.debug("OpenTargets: resolved %r -> %s", gene_name, target_id)
        return target_id

    # Search may return non-target IDs; try remaining hits
    for hit in hits[1:]:
        hit_id = hit.get("id", "")
        if hit_id.startswith("ENSG"):
            logger.debug("OpenTargets: resolved %r -> %s", gene_name, hit_id)
            return hit_id

    logger.debug(
        "OpenTargets: search for %r returned no ENSG IDs (first hit: %s)",
        gene_name,
        hits[0].get("id", ""),
    )
    return None


def _graphql_request(query: str, variables: dict[str, Any]) -> dict:
    """Execute a GraphQL POST request against the Open Targets API.

    Args:
        query: GraphQL query string.
        variables: Variables dict for the query.

    Returns:
        Parsed JSON response as a dict.

    Raises:
        Exception: On HTTP or parsing errors (after retries).
    """
    payload = {"query": query, "variables": variables}
    resp = post_with_retries(
        OPENTARGETS_GRAPHQL,
        json=payload,
        headers=_HEADERS,
    )
    return resp.json()
