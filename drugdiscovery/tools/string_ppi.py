"""STRING PPI network effect scoring.

Queries the STRING database API to build a local protein-protein
interaction network around the drug target, then propagates
perturbation effects through the network to estimate how many
disease-relevant genes are affected.

Tier 3 of the perturbation biology module (M4.6).
"""

from __future__ import annotations

import logging
from typing import Optional

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

STRING_API_BASE = "https://string-db.org/api"
SPECIES_HUMAN = 9606
MIN_INTERACTION_SCORE = 700  # STRING combined score threshold (0-1000)


# ---------------------------------------------------------------------------
# STRING API helpers
# ---------------------------------------------------------------------------

def _string_request(
    endpoint: str,
    params: dict,
    base_url: str = STRING_API_BASE,
    timeout: int = 30,
) -> Optional[list[dict]]:
    """Make a request to the STRING API.

    Returns parsed JSON list or None on failure.
    """
    import requests

    url = f"{base_url.rstrip('/')}/json/{endpoint}"

    try:
        resp = requests.get(url, params=params, timeout=timeout)
        resp.raise_for_status()
        data = resp.json()
        if isinstance(data, list):
            return data
        return [data] if isinstance(data, dict) else None
    except Exception as exc:
        logger.debug("[STRING] API request failed (%s): %s", endpoint, exc)
        return None


def resolve_string_id(
    gene_name: str,
    species: int = SPECIES_HUMAN,
    base_url: str = STRING_API_BASE,
) -> str:
    """Resolve a gene name to a STRING protein identifier.

    Returns the STRING ID (e.g. '9606.ENSP00000275493') or empty string.
    """
    result = _string_request(
        "get_string_ids",
        {
            "identifiers": gene_name,
            "species": species,
            "limit": 1,
            "echo_query": 1,
        },
        base_url=base_url,
    )
    if result and len(result) > 0:
        string_id = result[0].get("stringId", "")
        logger.debug("[STRING] Resolved %s → %s", gene_name, string_id)
        return string_id
    return ""


# ---------------------------------------------------------------------------
# Network retrieval
# ---------------------------------------------------------------------------

def get_interaction_partners(
    gene_name: str,
    species: int = SPECIES_HUMAN,
    limit: int = 50,
    min_score: int = MIN_INTERACTION_SCORE,
    base_url: str = STRING_API_BASE,
) -> list[dict]:
    """Get protein interaction partners from STRING.

    Returns list of dicts with:
      - partner_gene: str (gene symbol)
      - partner_string_id: str
      - combined_score: int (0-1000)
      - experimental_score: int
      - database_score: int
    """
    result = _string_request(
        "interaction_partners",
        {
            "identifiers": gene_name,
            "species": species,
            "limit": limit,
            "required_score": min_score,
        },
        base_url=base_url,
    )
    if not result:
        return []

    partners = []
    for item in result:
        # STRING returns preferredName_A/B and stringId_A/B
        # The query gene is usually in position A
        partner_name = item.get("preferredName_B", item.get("preferredName_A", ""))
        query_name = item.get("preferredName_A", "")

        # Make sure we pick the actual partner (not the query gene itself)
        if partner_name.upper() == gene_name.upper():
            partner_name = query_name

        partners.append({
            "partner_gene": partner_name,
            "partner_string_id": item.get("stringId_B", ""),
            "combined_score": int(item.get("score", 0) * 1000) if isinstance(item.get("score"), float) else int(item.get("combined_score", item.get("score", 0))),
            "experimental_score": int(item.get("escore", 0) * 1000) if isinstance(item.get("escore"), float) else int(item.get("experimental", 0)),
            "database_score": int(item.get("dscore", 0) * 1000) if isinstance(item.get("dscore"), float) else int(item.get("database", 0)),
        })

    # Deduplicate by partner name
    seen = set()
    unique = []
    for p in partners:
        name = p["partner_gene"].upper()
        if name and name not in seen and name != gene_name.upper():
            seen.add(name)
            unique.append(p)

    logger.info("[STRING] Found %d interaction partners for %s", len(unique), gene_name)
    return unique


def get_enrichment(
    gene_list: list[str],
    species: int = SPECIES_HUMAN,
    base_url: str = STRING_API_BASE,
) -> list[dict]:
    """Get functional enrichment for a gene list from STRING.

    Returns list of enriched terms with category, term, p-value.
    """
    if not gene_list:
        return []

    result = _string_request(
        "enrichment",
        {
            "identifiers": "%0d".join(gene_list[:200]),  # STRING limit
            "species": species,
        },
        base_url=base_url,
    )
    if not result:
        return []

    enrichments = []
    for item in result:
        enrichments.append({
            "category": item.get("category", ""),
            "term": item.get("term", ""),
            "description": item.get("description", ""),
            "p_value": float(item.get("p_value", 1.0)),
            "gene_count": int(item.get("number_of_genes", 0)),
            "genes": (
                item.get("inputGenes", [])
                if isinstance(item.get("inputGenes"), list)
                else item.get("inputGenes", "").split(",")
                if item.get("inputGenes")
                else []
            ),
        })

    # Sort by p-value
    enrichments.sort(key=lambda x: x["p_value"])
    return enrichments[:20]  # Top 20 enriched terms


# ---------------------------------------------------------------------------
# Network propagation scoring
# ---------------------------------------------------------------------------

def compute_network_effect(
    target_gene: str,
    disease_genes: list[str] | None = None,
    mode: str = "antagonist",
    depth: int = 2,
    base_url: str = STRING_API_BASE,
) -> dict:
    """Compute network propagation effect score.

    Models the downstream effect of perturbing the target gene by
    propagating through the STRING PPI network.

    Args:
        target_gene: Gene symbol of the drug target
        disease_genes: List of disease-relevant gene symbols
        mode: "agonist" or "antagonist"
        depth: Number of hops to propagate (1 or 2)
        base_url: STRING API base URL

    Returns dict with:
        - network_effect_score: float (0-1)
        - affected_genes: list[str]
        - affected_disease_genes: list[str]
        - total_network_genes: int
        - pathway_enrichment: list[dict]
    """
    disease_gene_set = set(g.upper() for g in (disease_genes or []))

    result = {
        "network_effect_score": 0.5,
        "affected_genes": [],
        "affected_disease_genes": [],
        "total_network_genes": 0,
        "pathway_enrichment": [],
    }

    # --- Depth 1: Direct interactors ---
    partners_d1 = get_interaction_partners(
        target_gene,
        limit=50,
        min_score=MIN_INTERACTION_SCORE,
        base_url=base_url,
    )

    if not partners_d1:
        logger.warning("[STRING] No interaction partners for %s; network_effect=0.5", target_gene)
        return result

    d1_genes = {p["partner_gene"].upper() for p in partners_d1}
    all_network_genes = d1_genes.copy()

    # --- Depth 2: Second-order interactors (optional) ---
    if depth >= 2:
        d2_genes: set[str] = set()
        # Only expand top 10 strongest interactors to limit API calls
        top_partners = sorted(partners_d1, key=lambda x: -x["combined_score"])[:10]
        for partner in top_partners:
            pname = partner["partner_gene"]
            if not pname:
                continue
            d2_partners = get_interaction_partners(
                pname,
                limit=20,
                min_score=MIN_INTERACTION_SCORE,
                base_url=base_url,
            )
            for p2 in d2_partners:
                g = p2["partner_gene"].upper()
                if g != target_gene.upper() and g not in d1_genes:
                    d2_genes.add(g)

        all_network_genes |= d2_genes

    result["affected_genes"] = sorted(all_network_genes)
    result["total_network_genes"] = len(all_network_genes)

    # --- Score: disease gene overlap ---
    if disease_gene_set:
        affected_disease = all_network_genes & disease_gene_set
        result["affected_disease_genes"] = sorted(affected_disease)

        # Score = fraction of disease genes affected by target perturbation
        if len(disease_gene_set) > 0:
            coverage = len(affected_disease) / len(disease_gene_set)
        else:
            coverage = 0.0

        # Boost score slightly if direct interactors overlap (stronger signal)
        d1_disease_overlap = d1_genes & disease_gene_set
        if d1_disease_overlap:
            d1_bonus = min(0.2, 0.05 * len(d1_disease_overlap))
            coverage = min(1.0, coverage + d1_bonus)

        result["network_effect_score"] = round(max(0.0, min(1.0, coverage)), 4)
    else:
        # No disease genes provided — use network size as proxy
        # More connections → more potential to affect biology
        # Normalize: 50+ interactors = score 0.8, fewer = proportional
        size_score = min(1.0, len(all_network_genes) / 60.0) * 0.8
        result["network_effect_score"] = round(max(0.1, size_score), 4)

    # --- Pathway enrichment ---
    if all_network_genes:
        enrichment = get_enrichment(
            list(all_network_genes)[:100],
            base_url=base_url,
        )
        result["pathway_enrichment"] = enrichment

    logger.info(
        "[STRING] Network effect for %s: score=%.4f, genes=%d, disease_overlap=%d",
        target_gene,
        result["network_effect_score"],
        result["total_network_genes"],
        len(result["affected_disease_genes"]),
    )

    return result


# ---------------------------------------------------------------------------
# Disease gene extraction helpers
# ---------------------------------------------------------------------------

def extract_disease_genes_from_go(go_terms: list[str]) -> list[str]:
    """Extract disease-relevant genes from GO biological process terms.

    Uses STRING's functional enrichment in reverse: given GO terms
    associated with the target, find genes in those pathways.

    This is a simplified approach; in production, one would use
    a dedicated disease gene database (DisGeNET, OMIM, etc.).
    """
    # For now, return empty — disease genes are better derived from
    # the target's interactome in the perturbation module itself
    return []


def get_target_pathway_genes(
    target_gene: str,
    base_url: str = STRING_API_BASE,
) -> list[str]:
    """Get all genes in pathways containing the target gene.

    Uses STRING enrichment to find KEGG/Reactome pathways, then
    returns the union of genes in those pathways.
    """
    enrichment = get_enrichment(
        [target_gene],
        base_url=base_url,
    )

    pathway_genes: set[str] = set()
    for term in enrichment:
        if term["category"] in ("KEGG", "Reactome", "WikiPathways"):
            pathway_genes.update(g.strip() for g in term.get("genes", []))

    logger.debug(
        "[STRING] Found %d pathway genes for %s across %d terms",
        len(pathway_genes), target_gene, len(enrichment),
    )
    return sorted(pathway_genes)
