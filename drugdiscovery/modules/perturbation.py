"""M4.6 – Perturbation Biology Scoring module.

Predicts how each drug candidate will perturb disease-relevant gene
expression at the cellular level.  Combines three tiers:

  T1: L1000/CMAP connectivity scoring (weight 0.50)
  T2: Tanimoto-kNN perturbation transfer (weight 0.20)
  T3: STRING PPI network effect (weight 0.30)

Produces a perturbation_score (0-1) for each candidate, where higher
indicates stronger evidence that the candidate will produce the desired
cellular effect (disease signature reversal / on-target network effect).
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

from drugdiscovery.types import Candidate, PipelineConfig, TargetProfile
from drugdiscovery.utils.io import write_csv

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Default tier weights
# ---------------------------------------------------------------------------

DEFAULT_TIER_WEIGHTS = {
    "cmap_connectivity": 0.50,
    "knn_transfer": 0.20,
    "network_effect": 0.30,
}


# ---------------------------------------------------------------------------
# Disease gene derivation
# ---------------------------------------------------------------------------

def _derive_disease_genes(
    target_profile: TargetProfile,
    string_base_url: str = "https://string-db.org/api",
) -> tuple[list[str], list[str]]:
    """Derive disease-relevant up/down gene sets from the target profile.

    Strategy:
      1. Get first-order interactors from STRING PPI
      2. Genes in the same GO biological processes as the target
      3. Anti-target genes as potential off-target markers

    Returns (disease_up_genes, disease_down_genes).
    For an antagonist target, genes upregulated in disease are those
    that the target activates (interactors).
    """
    disease_up: list[str] = []
    disease_down: list[str] = []

    try:
        from drugdiscovery.tools.string_ppi import get_interaction_partners

        partners = get_interaction_partners(
            target_profile.gene_name,
            limit=30,
            base_url=string_base_url,
        )
        # Target interactors are assumed upregulated in the disease context
        # (the target drives their expression)
        disease_up = [p["partner_gene"] for p in partners if p["partner_gene"]]

    except Exception as exc:
        logger.warning("[M4.6] Failed to get STRING interactors: %s", exc)

    # Add anti-target genes to disease_down (genes that should NOT be affected)
    for anti in target_profile.anti_targets:
        gene = anti.get("gene_name", "")
        if gene:
            disease_down.append(gene)

    # Also derive from GO terms if available
    if target_profile.go_terms:
        # Use GO terms as keywords — these genes are functionally related
        # to the target and likely affected in the disease
        for term in target_profile.go_terms[:5]:
            # Extract gene-like tokens (this is heuristic)
            if ":" in term:
                continue  # Skip GO IDs like "GO:0006412"
            # Some GO terms contain gene names
            pass  # Placeholder for future GO-based gene extraction

    logger.info(
        "[M4.6] Derived disease genes: %d up, %d down",
        len(disease_up), len(disease_down),
    )
    return disease_up, disease_down


# ---------------------------------------------------------------------------
# Per-candidate perturbation scoring
# ---------------------------------------------------------------------------

def _score_candidate_perturbation(
    candidate: Candidate,
    target_profile: TargetProfile,
    disease_up_genes: list[str],
    disease_down_genes: list[str],
    network_result: dict,
    tier_weights: dict,
    cmap_api_key: str = "",
    cmap_base_url: str = "https://api.clue.io/api",
    knn_k: int = 5,
) -> Candidate:
    """Compute perturbation score for a single candidate.

    Combines three tiers into a weighted aggregate.
    """
    from drugdiscovery.tools.cmap_client import query_cmap_perturbation

    w = {**DEFAULT_TIER_WEIGHTS, **tier_weights}
    total_weight = sum(w.values())

    # --- Tier 1 + 2: CMAP connectivity + kNN transfer ---
    cmap_result = query_cmap_perturbation(
        smiles=candidate.smiles,
        compound_name=candidate.source if candidate.source not in ("pepmlm", "de_novo", "diffpepbuilder") else "",
        gene_target=target_profile.gene_name,
        disease_up_genes=disease_up_genes,
        disease_down_genes=disease_down_genes,
        api_key=cmap_api_key,
        base_url=cmap_base_url,
        knn_k=knn_k,
    )

    cmap_connectivity = cmap_result.get("connectivity_score", 0.5)
    matched_compound = cmap_result.get("matched_compound", "")
    match_source = cmap_result.get("match_source", "none")

    # For kNN transfer tier: if match came from kNN, use a discounted score
    if match_source == "knn_transfer":
        knn_score = cmap_connectivity * 0.8  # Discount transferred scores
        cmap_score = 0.5  # No direct CMAP data
    elif match_source == "clue_api":
        cmap_score = cmap_connectivity
        knn_score = cmap_connectivity  # Direct match is best kNN too
    else:
        cmap_score = 0.5
        knn_score = 0.5

    # --- Tier 3: Network effect ---
    network_score = network_result.get("network_effect_score", 0.5)

    # --- Peptide-specific adjustments ---
    if candidate.modality == "peptide":
        # Peptides are less likely to be in L1000 (small molecule focused)
        # Rely more on network effect for peptides
        w["cmap_connectivity"] = w.get("cmap_connectivity", 0.5) * 0.5
        w["network_effect"] = w.get("network_effect", 0.3) * 1.5
        total_weight = sum(w.values())

    # --- Weighted aggregate ---
    if total_weight < 1e-9:
        perturbation_score = 0.5
    else:
        perturbation_score = (
            w.get("cmap_connectivity", 0.5) * cmap_score
            + w.get("knn_transfer", 0.2) * knn_score
            + w.get("network_effect", 0.3) * network_score
        ) / total_weight

    perturbation_score = round(max(0.0, min(1.0, perturbation_score)), 4)

    # --- Assign to candidate ---
    candidate.perturbation_score = perturbation_score
    candidate.cmap_connectivity = round(cmap_connectivity, 4)
    candidate.cmap_compound_match = matched_compound
    candidate.network_effect_score = round(network_score, 4)
    candidate.disease_signature_reversal = round(cmap_score, 4)

    logger.debug(
        "[M4.6] %s perturbation=%.4f (cmap=%.4f knn=%.4f net=%.4f match=%s)",
        candidate.candidate_id,
        perturbation_score,
        cmap_score,
        knn_score,
        network_score,
        matched_compound or "none",
    )

    return candidate


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def predict_perturbation(
    cfg: PipelineConfig,
    candidates: list[Candidate],
    target_profile: TargetProfile,
    output_dir: Path,
) -> list[Candidate]:
    """M4.6 — Predict cellular perturbation effects for all candidates.

    Steps:
      1. Derive disease gene sets from target profile + STRING PPI
      2. Compute network effect (once per target, shared across candidates)
      3. For each candidate, compute perturbation score via 3 tiers
      4. Save perturbation results to output_dir

    Returns candidates with perturbation_score populated.
    """
    if not candidates:
        logger.warning("[M4.6] No candidates to score.")
        return candidates

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Load perturbation config ---
    pert_cfg = getattr(cfg, "perturbation_settings", {})
    if not pert_cfg and hasattr(cfg, "config_file") and cfg.config_file:
        # Try loading from YAML config
        try:
            import yaml
            with open(cfg.config_file) as f:
                all_cfg = yaml.safe_load(f) or {}
            pert_cfg = all_cfg.get("perturbation", {})
        except Exception:
            pass

    cmap_api_key = pert_cfg.get("cmap_api_key", "")
    cmap_base_url = pert_cfg.get("cmap_api_url", "https://api.clue.io/api")
    string_base_url = pert_cfg.get("string_api_url", "https://string-db.org/api")
    network_depth = pert_cfg.get("network_depth", 2)
    knn_k = pert_cfg.get("knn_neighbors", 5)
    tier_weights = pert_cfg.get("tier_weights", DEFAULT_TIER_WEIGHTS)

    logger.info("[M4.6] Perturbation scoring for %d candidates", len(candidates))
    logger.info("[M4.6] Target: %s, Tier weights: %s", target_profile.gene_name, tier_weights)

    # --- Step 1: Derive disease gene sets ---
    disease_up, disease_down = _derive_disease_genes(
        target_profile,
        string_base_url=string_base_url,
    )

    # --- Step 2: Network effect (computed once for the target) ---
    network_result = {"network_effect_score": 0.5}
    try:
        from drugdiscovery.tools.string_ppi import compute_network_effect

        network_result = compute_network_effect(
            target_gene=target_profile.gene_name,
            disease_genes=disease_up + disease_down,
            mode=cfg.mode.value,
            depth=network_depth,
            base_url=string_base_url,
        )
        logger.info(
            "[M4.6] Network effect: score=%.4f, %d genes, %d disease overlap",
            network_result["network_effect_score"],
            network_result.get("total_network_genes", 0),
            len(network_result.get("affected_disease_genes", [])),
        )
    except Exception as exc:
        logger.warning("[M4.6] Network effect computation failed: %s", exc, exc_info=True)

    # --- Step 3: Per-candidate scoring ---
    for i, cand in enumerate(candidates):
        try:
            _score_candidate_perturbation(
                cand,
                target_profile,
                disease_up,
                disease_down,
                network_result,
                tier_weights,
                cmap_api_key=cmap_api_key,
                cmap_base_url=cmap_base_url,
                knn_k=knn_k,
            )
        except Exception as exc:
            logger.warning(
                "[M4.6] Perturbation scoring failed for %s: %s",
                cand.candidate_id, exc,
            )
            cand.perturbation_score = 0.5

        if (i + 1) % 50 == 0:
            logger.info("[M4.6] Scored %d/%d candidates", i + 1, len(candidates))

    # --- Step 4: Save results ---
    pert_rows = []
    for cand in candidates:
        pert_rows.append({
            "candidate_id": cand.candidate_id,
            "modality": cand.modality,
            "perturbation_score": round(cand.perturbation_score, 4),
            "cmap_connectivity": round(cand.cmap_connectivity, 4),
            "cmap_compound_match": cand.cmap_compound_match,
            "network_effect_score": round(cand.network_effect_score, 4),
            "disease_signature_reversal": round(cand.disease_signature_reversal, 4),
        })

    write_csv(output_dir / "perturbation_scores.csv", pert_rows)
    logger.info("[M4.6] Saved perturbation_scores.csv (%d rows)", len(pert_rows))

    # Save network effect details
    try:
        with open(output_dir / "network_effect.json", "w") as f:
            # Convert sets to lists for JSON serialization
            serializable = {
                k: (list(v) if isinstance(v, set) else v)
                for k, v in network_result.items()
            }
            json.dump(serializable, f, indent=2, default=str)
        logger.info("[M4.6] Saved network_effect.json")
    except Exception as exc:
        logger.warning("[M4.6] Failed to save network_effect.json: %s", exc)

    # Save disease genes used
    try:
        with open(output_dir / "disease_genes.json", "w") as f:
            json.dump({
                "target": target_profile.gene_name,
                "disease_up_genes": disease_up,
                "disease_down_genes": disease_down,
            }, f, indent=2)
    except Exception:
        pass

    # Summary statistics
    scores = [c.perturbation_score for c in candidates]
    avg_score = sum(scores) / len(scores) if scores else 0.0
    logger.info(
        "[M4.6] Perturbation scoring complete. "
        "Mean=%.4f, Min=%.4f, Max=%.4f",
        avg_score, min(scores), max(scores),
    )

    return candidates
