"""Pipeline orchestrator: wires all modules together."""

from __future__ import annotations

import logging
import time
from pathlib import Path

from drugdiscovery.config import (
    make_run_dir,
    resolve_device,
    save_config_snapshot,
    setup_logging,
)
from drugdiscovery.types import Candidate, Modality, PipelineConfig

logger = logging.getLogger(__name__)


def run_pipeline(cfg: PipelineConfig) -> Path:
    """Execute the full drug discovery pipeline.

    Returns the path to the run output directory.
    """
    # --- Setup ---
    device = resolve_device(cfg.device)
    cfg.device = device
    run_dir = make_run_dir(cfg.output_dir, cfg.target, cfg.modality.value, cfg.mode.value)
    setup_logging(run_dir)
    save_config_snapshot(cfg, run_dir)

    logger.info("=" * 70)
    logger.info("DRUG DISCOVERY PIPELINE")
    logger.info("  Target:   %s", cfg.target)
    logger.info("  Modality: %s", cfg.modality.value)
    logger.info("  Mode:     %s", cfg.mode.value)
    logger.info("  Device:   %s", device)
    logger.info("  Output:   %s", run_dir)
    logger.info("=" * 70)

    t0 = time.time()

    # --- M1: Target Preparation ---
    logger.info("[M1] Target Preparation")
    from drugdiscovery.modules.target_prep import prepare_target

    target_profile = prepare_target(cfg, run_dir / "target")
    logger.info("[M1] Target resolved: %s (%s)", target_profile.gene_name, target_profile.uniprot_id)

    # --- M2 + M3: Library Screening + De Novo Design ---
    candidates: list[Candidate] = []

    # M2: Library Screening
    logger.info("[M2] Library Screening")
    try:
        from drugdiscovery.modules.library_screening import screen_libraries

        library_hits = screen_libraries(cfg, target_profile, run_dir / "library_hits")
        candidates.extend(library_hits)
        logger.info("[M2] Library hits: %d candidates", len(library_hits))
    except Exception as exc:
        logger.warning("[M2] Library screening failed (continuing): %s", exc, exc_info=True)

    # M3: De Novo Design
    logger.info("[M3] De Novo Design")
    if cfg.modality == Modality.PEPTIDE:
        from drugdiscovery.modules.de_novo_peptide import design_peptides

        de_novo = design_peptides(cfg, target_profile, run_dir / "de_novo")
    else:
        from drugdiscovery.modules.de_novo_sm import design_small_molecules

        de_novo = design_small_molecules(cfg, target_profile, run_dir / "de_novo")
    candidates.extend(de_novo)
    logger.info("[M3] De novo candidates: %d", len(de_novo))

    if not candidates:
        logger.error("No candidates generated. Pipeline cannot continue.")
        return run_dir

    # --- M4: Structure Prediction ---
    logger.info("[M4] Structure Prediction")
    try:
        from drugdiscovery.modules.structure_pred import predict_structures

        candidates = predict_structures(cfg, candidates, target_profile, run_dir)
        logger.info("[M4] Structure prediction complete for %d candidates", len(candidates))
    except Exception as exc:
        logger.warning("[M4] Structure prediction skipped: %s", exc, exc_info=True)

    # --- M4.5: Molecular Docking ---
    logger.info("[M4.5] Molecular Docking")
    try:
        from drugdiscovery.modules.docking import dock_candidates

        candidates = dock_candidates(cfg, candidates, target_profile, run_dir)
        logger.info("[M4.5] Docking complete for %d candidates", len(candidates))
    except Exception as exc:
        logger.warning("[M4.5] Docking skipped: %s", exc, exc_info=True)

    # --- M7: ADMET Prediction (before scoring so admet_score enters composite) ---
    logger.info("[M7] ADMET Prediction")
    from drugdiscovery.modules.admet import predict_admet

    candidates = predict_admet(cfg, candidates, run_dir / "admet")
    logger.info("[M7] ADMET profiles computed for %d candidates", len(candidates))

    # --- M5: Binding & Scoring ---
    logger.info("[M5] Binding & Scoring")
    from drugdiscovery.modules.scoring import score_candidates

    candidates = score_candidates(cfg, candidates, target_profile, run_dir / "scoring")
    logger.info("[M5] Scored %d candidates", len(candidates))

    # Filter to top N before modifications
    candidates.sort(key=lambda c: c.composite_score, reverse=True)
    top_candidates = candidates[: cfg.top_n]
    logger.info("[M5] Top %d candidates selected", len(top_candidates))

    # --- M6: Modifications ---
    logger.info("[M6] Modifications")
    from drugdiscovery.modules.modifications import modify_candidates

    modified = modify_candidates(cfg, top_candidates, run_dir / "modifications")
    all_candidates = top_candidates + modified
    logger.info("[M6] Generated %d modified variants", len(modified))

    # Re-dock modified candidates so binding_score reflects modifications
    if modified:
        logger.info("[M4.5] Re-docking %d modified candidates", len(modified))
        try:
            from drugdiscovery.modules.docking import dock_candidates  # noqa: PLC0415

            modified = dock_candidates(cfg, modified, target_profile, run_dir)
            logger.info("[M4.5] Re-docking complete for %d modified candidates", len(modified))
            all_candidates = top_candidates + modified
        except Exception as exc:
            logger.warning("[M4.5] Re-docking skipped: %s", exc, exc_info=True)

    # Re-run ADMET for modified candidates
    logger.info("[M7] ADMET for modified candidates")
    all_candidates = predict_admet(cfg, all_candidates, run_dir / "admet")

    # Re-score all candidates with updated binding + ADMET
    all_candidates = score_candidates(cfg, all_candidates, target_profile, run_dir / "scoring")

    # --- M8: Delivery System ---
    logger.info("[M8] Delivery System Recommendation")
    from drugdiscovery.modules.delivery import recommend_delivery

    all_candidates = recommend_delivery(cfg, all_candidates, target_profile, run_dir / "delivery")
    logger.info("[M8] Delivery recommendations assigned")

    # --- Final ranking ---
    all_candidates.sort(key=lambda c: c.composite_score, reverse=True)
    for i, c in enumerate(all_candidates, 1):
        c.rank = i

    # --- M9: Final Report ---
    logger.info("[M9] Final Report")
    from drugdiscovery.modules.reporting import generate_report

    generate_report(cfg, all_candidates, target_profile, run_dir / "report")

    elapsed = time.time() - t0
    logger.info("=" * 70)
    logger.info("PIPELINE COMPLETE")
    logger.info("  Total candidates: %d", len(all_candidates))
    logger.info("  Top candidate: %s (score=%.4f)", all_candidates[0].candidate_id, all_candidates[0].composite_score)
    logger.info("  Elapsed: %.1f seconds", elapsed)
    logger.info("  Results: %s", run_dir)
    logger.info("=" * 70)

    return run_dir
