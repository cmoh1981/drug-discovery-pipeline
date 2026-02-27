"""M2: De Novo Peptide Design.

Generates peptide candidates using two complementary approaches:
  1. PepMLM  – sequence-based masked-language-model generation; always runs
               and degrades gracefully to CPU when no GPU is available.
  2. DiffPepBuilder – structure-based diffusion design; runs only when a
                      structure PDB is available in the TargetProfile.

Results from both tools are merged, deduplicated by sequence, converted to
Candidate objects, and saved as candidates.csv and candidates.fasta.
"""

from __future__ import annotations

import logging
import math
import uuid
from pathlib import Path

from drugdiscovery.types import Candidate, PipelineConfig, TargetProfile
from drugdiscovery.utils.io import write_csv, write_fasta
from drugdiscovery.tools import pepmlm as _pepmlm_mod
from drugdiscovery.tools import diffpepbuilder as _dpb_mod

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _make_candidate_id(prefix: str) -> str:
    """Generate a short unique candidate ID with a readable prefix."""
    short = uuid.uuid4().hex[:8].upper()
    return f"{prefix}_{short}"


def _is_valid_sequence(seq: str) -> bool:
    """Return True if *seq* is a non-empty string of standard amino acids."""
    if not seq:
        return False
    return all(aa in _pepmlm_mod.AA_TOKENS for aa in seq.upper())


def _perplexity_to_float(value: object) -> float:
    """Coerce a perplexity value to float; return NaN on failure."""
    if isinstance(value, float):
        return value
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return float("nan")


# ---------------------------------------------------------------------------
# PepMLM runner
# ---------------------------------------------------------------------------

def _run_pepmlm(
    cfg: PipelineConfig,
    target_profile: TargetProfile,
) -> list[dict]:
    """Run PepMLM across all configured lengths and return raw result dicts.

    On any failure the function logs a warning and returns an empty list so
    that the orchestrator can continue with DiffPepBuilder results only.

    Returns:
        List of dicts with keys: ``sequence``, ``perplexity``, ``length``.
    """
    settings = cfg.pepmlm_settings
    model_name: str = settings.get("model_name", _pepmlm_mod.PRIMARY_MODEL)
    lengths: list[int] = settings.get("lengths", [8, 10, 12, 15, 20, 25, 30])
    num_per_length: int = settings.get("num_per_length", 50)
    top_k: int = settings.get("top_k", 3)

    # Resolve device – honour cfg.device; fall back to CPU so PepMLM always runs.
    from drugdiscovery.utils.compute import detect_device  # noqa: PLC0415
    device = detect_device(cfg.device)

    target_sequence = target_profile.sequence
    if not target_sequence:
        logger.warning("[M2] TargetProfile has no sequence; skipping PepMLM.")
        return []

    logger.info(
        "[M2] PepMLM | model=%s  device=%s  lengths=%s  num_per_length=%d  top_k=%d",
        model_name,
        device,
        lengths,
        num_per_length,
        top_k,
    )

    try:
        model, tokenizer, device = _pepmlm_mod.load_pepmlm(model_name, device)
    except Exception as exc:  # noqa: BLE001
        logger.warning("[M2] PepMLM model load failed – skipping PepMLM: %s", exc)
        return []

    raw: list[dict] = []
    for length in lengths:
        logger.info("[M2] PepMLM | generating %d peptides of length %d …", num_per_length, length)
        try:
            batch = _pepmlm_mod.generate_peptides(
                model=model,
                tokenizer=tokenizer,
                device=device,
                target_sequence=target_sequence,
                peptide_length=length,
                num_peptides=num_per_length,
                top_k=top_k,
            )
        except Exception as exc:  # noqa: BLE001
            logger.warning(
                "[M2] PepMLM generation failed for length %d: %s", length, exc
            )
            continue

        for item in batch:
            seq = item.get("sequence", "")
            if _is_valid_sequence(seq):
                raw.append({
                    "sequence": seq.upper(),
                    "perplexity": item.get("perplexity", float("nan")),
                    "length": len(seq),
                })

    logger.info("[M2] PepMLM | produced %d valid sequences across %d lengths.", len(raw), len(lengths))
    return raw


# ---------------------------------------------------------------------------
# DiffPepBuilder runner
# ---------------------------------------------------------------------------

def _run_diffpepbuilder(
    cfg: PipelineConfig,
    target_profile: TargetProfile,
    output_dir: Path,
) -> list[dict]:
    """Run DiffPepBuilder (local or RunPod) and return raw result dicts.

    Skipped silently when no structure PDB is available.  On any other failure
    a warning is logged and an empty list is returned.

    Returns:
        List of dicts with keys: ``sequence``, ``pdb_path``, ``confidence``.
    """
    if not target_profile.structure_pdb_path:
        logger.info(
            "[M2] DiffPepBuilder | no structure PDB available for '%s'; skipping.",
            target_profile.gene_name,
        )
        return []

    settings = cfg.diffpepbuilder_settings
    num_samples: int = settings.get("num_samples", 100)
    pep_min: int = settings.get("peptide_length_min", 8)
    pep_max: int = settings.get("peptide_length_max", 25)

    # Collect binding site residue numbers from all known pockets
    binding_residues: list[int] = []
    for pocket in target_profile.binding_pockets:
        binding_residues.extend(pocket.residue_numbers)
    binding_residues = sorted(set(binding_residues))

    if not binding_residues:
        logger.warning(
            "[M2] DiffPepBuilder | no binding-site residues defined for '%s'; skipping.",
            target_profile.gene_name,
        )
        return []

    dpb_output = output_dir / "diffpepbuilder_output"
    dpb_output.mkdir(parents=True, exist_ok=True)

    logger.info(
        "[M2] DiffPepBuilder | pdb=%s  hotspots=%d residues  samples=%d  len=[%d,%d]  runpod=%s",
        target_profile.structure_pdb_path,
        len(binding_residues),
        num_samples,
        pep_min,
        pep_max,
        cfg.use_runpod,
    )

    try:
        if cfg.use_runpod:
            raw = _dpb_mod.run_diffpepbuilder_runpod(
                target_pdb=target_profile.structure_pdb_path,
                binding_site_residues=binding_residues,
                output_dir=str(dpb_output),
                num_samples=num_samples,
                peptide_length_min=pep_min,
                peptide_length_max=pep_max,
            )
        else:
            raw = _dpb_mod.run_diffpepbuilder_local(
                target_pdb=target_profile.structure_pdb_path,
                binding_site_residues=binding_residues,
                output_dir=str(dpb_output),
                num_samples=num_samples,
                peptide_length_min=pep_min,
                peptide_length_max=pep_max,
            )
    except Exception as exc:  # noqa: BLE001
        logger.warning(
            "[M2] DiffPepBuilder run failed – continuing without it: %s", exc
        )
        return []

    # Filter to valid amino acid sequences
    valid = [r for r in raw if _is_valid_sequence(r.get("sequence", ""))]
    logger.info(
        "[M2] DiffPepBuilder | %d/%d results passed sequence validation.",
        len(valid),
        len(raw),
    )
    return valid


# ---------------------------------------------------------------------------
# Candidate construction
# ---------------------------------------------------------------------------

def _pepmlm_to_candidate(item: dict) -> Candidate:
    """Convert a raw PepMLM result dict to a Candidate."""
    seq = item["sequence"].upper()
    ppl = _perplexity_to_float(item.get("perplexity"))
    return Candidate(
        candidate_id=_make_candidate_id("PEPMLM"),
        candidate_type="de_novo",
        modality="peptide",
        source="pepmlm",
        sequence=seq,
        perplexity=ppl if not math.isnan(ppl) else 0.0,
    )


def _dpb_to_candidate(item: dict) -> Candidate:
    """Convert a raw DiffPepBuilder result dict to a Candidate."""
    seq = item["sequence"].upper()
    try:
        conf = float(item.get("confidence", float("nan")))
    except (TypeError, ValueError):
        conf = float("nan")

    return Candidate(
        candidate_id=_make_candidate_id("DPB"),
        candidate_type="de_novo",
        modality="peptide",
        source="diffpepbuilder",
        sequence=seq,
        structure_confidence=conf if not math.isnan(conf) else 0.0,
        metadata={"pdb_path": item.get("pdb_path", "")},
    )


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def design_peptides(
    cfg: PipelineConfig,
    target_profile: TargetProfile,
    output_dir: Path,
) -> list[Candidate]:
    """Generate de novo peptide candidates using PepMLM and DiffPepBuilder.

    Workflow:
      1. Run PepMLM for sequence-based generation (always attempted; falls back
         to CPU if no GPU; returns empty list on any tool failure).
      2. Run DiffPepBuilder for structure-based generation (skipped when no
         structure PDB is present; returns empty list on tool failure).
      3. Merge results; deduplicate by exact sequence (first occurrence wins).
      4. Convert to Candidate objects.
      5. Save candidates.csv and candidates.fasta to *output_dir*.
      6. Return the Candidate list.

    Args:
        cfg:            Pipeline configuration.
        target_profile: Resolved target protein profile.
        output_dir:     Directory for all M2 artefacts.

    Returns:
        List of Candidate objects (may be empty if both tools fail).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(
        "[M2] De novo peptide design | target=%s  runpod=%s",
        target_profile.gene_name,
        cfg.use_runpod,
    )

    # --- Step 1: PepMLM ---
    pepmlm_results: list[dict] = []
    try:
        pepmlm_results = _run_pepmlm(cfg, target_profile)
    except Exception as exc:  # noqa: BLE001
        logger.warning("[M2] PepMLM raised an unexpected error – skipping: %s", exc)

    if not pepmlm_results:
        logger.warning(
            "[M2] PepMLM produced no candidates for '%s'.", target_profile.gene_name
        )

    # --- Step 2: DiffPepBuilder ---
    dpb_results: list[dict] = []
    try:
        dpb_results = _run_diffpepbuilder(cfg, target_profile, output_dir)
    except Exception as exc:  # noqa: BLE001
        logger.warning(
            "[M2] DiffPepBuilder raised an unexpected error – skipping: %s", exc
        )

    # --- Step 3: Merge and deduplicate by sequence ---
    seen_sequences: set[str] = set()
    candidates: list[Candidate] = []

    for item in pepmlm_results:
        seq = item["sequence"].upper()
        if seq in seen_sequences:
            continue
        seen_sequences.add(seq)
        candidates.append(_pepmlm_to_candidate(item))

    for item in dpb_results:
        seq = item["sequence"].upper()
        if seq in seen_sequences:
            logger.debug("[M2] Duplicate sequence from DiffPepBuilder skipped: %s", seq)
            continue
        seen_sequences.add(seq)
        candidates.append(_dpb_to_candidate(item))

    logger.info(
        "[M2] Merged candidates: %d PepMLM + %d DiffPepBuilder = %d unique.",
        len(pepmlm_results),
        len(dpb_results),
        len(candidates),
    )

    if not candidates:
        logger.warning(
            "[M2] No de novo candidates generated for '%s'. "
            "Check tool availability, model access, and target structure.",
            target_profile.gene_name,
        )
        return candidates

    # --- Step 4: Save candidates.csv ---
    csv_rows = [c.to_dict() for c in candidates]
    write_csv(output_dir / "candidates.csv", csv_rows)
    logger.info("[M2] Saved candidates.csv (%d rows).", len(csv_rows))

    # --- Step 5: Save candidates.fasta ---
    fasta_entries: list[tuple[str, str]] = []
    for c in candidates:
        ppl_str = f"{c.perplexity:.4f}" if c.perplexity else "NA"
        conf_str = f"{c.structure_confidence:.4f}" if c.structure_confidence else "NA"
        header = (
            f"{c.candidate_id} "
            f"source={c.source} "
            f"perplexity={ppl_str} "
            f"structure_confidence={conf_str}"
        )
        fasta_entries.append((header, c.sequence))
    write_fasta(output_dir / "candidates.fasta", fasta_entries)
    logger.info("[M2] Saved candidates.fasta (%d sequences).", len(fasta_entries))

    logger.info(
        "[M2] De novo peptide design complete: %d candidates "
        "(%d from PepMLM, %d from DiffPepBuilder).",
        len(candidates),
        sum(1 for c in candidates if c.source == "pepmlm"),
        sum(1 for c in candidates if c.source == "diffpepbuilder"),
    )

    return candidates
