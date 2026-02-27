"""PepMLM wrapper for sequence-based peptide binder generation.

Refactored from YARS2 pipeline: scripts/run_pepmlm.py

Model: TianlaiChen/PepMLM-650M (fallback: ChatterjeeLab/PepMLM-650M)

All heavy ML imports (torch, transformers) are lazy so the module can be
imported without them installed.  Missing-library errors surface only when
the generation functions are actually called.
"""

from __future__ import annotations

import logging
import math
from typing import Any

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PRIMARY_MODEL = "TianlaiChen/PepMLM-650M"
FALLBACK_MODEL = "ChatterjeeLab/PepMLM-650M"

# Standard 20 amino acids – generation and scoring are restricted to these.
AA_TOKENS: set[str] = set("ACDEFGHIKLMNPQRSTVWY")


# ---------------------------------------------------------------------------
# Lazy imports
# ---------------------------------------------------------------------------

def _import_torch() -> Any:
    """Lazily import torch; raises ImportError with a clear message if absent."""
    try:
        import torch  # noqa: PLC0415
        return torch
    except ImportError as exc:
        raise ImportError(
            "torch is required for PepMLM generation. "
            "Install it with: pip install torch"
        ) from exc


def _import_transformers() -> tuple[Any, Any]:
    """Lazily import AutoTokenizer and AutoModelForMaskedLM."""
    try:
        from transformers import AutoTokenizer, AutoModelForMaskedLM  # noqa: PLC0415
        return AutoTokenizer, AutoModelForMaskedLM
    except ImportError as exc:
        raise ImportError(
            "transformers is required for PepMLM generation. "
            "Install it with: pip install transformers"
        ) from exc


# ---------------------------------------------------------------------------
# Model loading
# ---------------------------------------------------------------------------

def load_pepmlm(
    model_name: str = PRIMARY_MODEL,
    device: str = "cpu",
) -> tuple[Any, Any, str]:
    """Load PepMLM model and tokenizer.

    Tries *model_name* first; if that fails and *model_name* is the primary
    HuggingFace checkpoint, automatically retries with the fallback checkpoint.

    Args:
        model_name: HuggingFace model identifier.
        device:     Torch device string ("cpu", "cuda", "mps").

    Returns:
        (model, tokenizer, device_str) – the device_str reflects the actual
        device used (may differ from *device* if the requested device is
        unavailable).

    Raises:
        RuntimeError: If no model variant can be loaded.
        ImportError:  If torch or transformers are not installed.
    """
    torch = _import_torch()
    AutoTokenizer, AutoModelForMaskedLM = _import_transformers()

    # Build fallback chain
    names_to_try: list[str] = [model_name]
    if model_name == PRIMARY_MODEL:
        names_to_try.append(FALLBACK_MODEL)

    last_exc: Exception | None = None
    for name in names_to_try:
        try:
            logger.info("[PepMLM] Loading tokenizer from '%s' …", name)
            tokenizer = AutoTokenizer.from_pretrained(name)
            logger.info("[PepMLM] Loading model from '%s' …", name)
            model = AutoModelForMaskedLM.from_pretrained(name)
            model = model.to(device)
            model.eval()
            logger.info("[PepMLM] Model ready on device '%s'.", device)
            return model, tokenizer, device
        except Exception as exc:  # noqa: BLE001
            logger.warning("[PepMLM] Could not load '%s': %s", name, exc)
            last_exc = exc

    raise RuntimeError(
        f"[PepMLM] Could not load any model from {names_to_try}. "
        f"Last error: {last_exc}"
    )


# ---------------------------------------------------------------------------
# Amino-acid token ID pre-computation (module-level cache per tokenizer)
# ---------------------------------------------------------------------------

def _get_aa_token_ids(tokenizer: Any) -> list[int]:
    """Return the vocabulary token IDs that correspond to standard amino acids."""
    aa_token_ids: list[int] = []
    for aa in AA_TOKENS:
        ids = tokenizer.encode(aa, add_special_tokens=False)
        if ids:
            aa_token_ids.append(ids[0])
    return aa_token_ids


def _get_mask_token_id(tokenizer: Any) -> int:
    """Return the mask token ID, trying both .mask_token_id and <mask>."""
    mid = tokenizer.mask_token_id
    if mid is None:
        mid = tokenizer.convert_tokens_to_ids("<mask>")
    return mid


# ---------------------------------------------------------------------------
# Peptide generation
# ---------------------------------------------------------------------------

def generate_peptides(
    model: Any,
    tokenizer: Any,
    device: str,
    target_sequence: str,
    peptide_length: int = 15,
    num_peptides: int = 50,
    top_k: int = 3,
) -> list[dict]:
    """Generate peptide binder candidates using masked-language modelling.

    Input construction:
        ``[CLS] target_sequence [SEP] [MASK]*peptide_length [SEP]``

    For each of *num_peptides* independent samples:
      - Run a single forward pass to get logits at all MASK positions.
      - At each masked position, zero out non-amino-acid tokens, apply
        top-k filtering, and sample one amino acid from the softmax.

    Args:
        model:           Loaded PepMLM model (eval mode, on *device*).
        tokenizer:       Corresponding HuggingFace tokenizer.
        device:          Torch device string.
        target_sequence: Full amino acid sequence of the target protein.
        peptide_length:  Desired peptide length in residues.
        num_peptides:    Number of independent candidate sequences to generate.
        top_k:           Top-k for amino acid sampling at each masked position.

    Returns:
        List of dicts with keys:
          ``sequence``   – generated amino acid sequence (str)
          ``perplexity`` – pseudo-perplexity score (float); NaN if scoring
                           failed (individual scoring errors are non-fatal).
    """
    torch = _import_torch()

    aa_token_ids = _get_aa_token_ids(tokenizer)
    if not aa_token_ids:
        raise RuntimeError(
            "[PepMLM] Could not resolve any amino-acid token IDs from tokenizer vocabulary. "
            "Verify that the loaded model is an ESM-family protein language model."
        )

    mask_token_id = _get_mask_token_id(tokenizer)
    sep = tokenizer.eos_token if tokenizer.eos_token else " "
    mask_str = tokenizer.mask_token * peptide_length

    candidates: list[dict] = []

    with torch.no_grad():
        for sample_idx in range(num_peptides):
            concat_input = target_sequence + sep + mask_str
            inputs = tokenizer(
                concat_input,
                return_tensors="pt",
                truncation=True,
                max_length=tokenizer.model_max_length,
            )
            inputs = {k: v.to(device) for k, v in inputs.items()}
            input_ids = inputs["input_ids"]  # [1, seq_len]

            # Locate masked positions
            mask_positions = (
                input_ids[0] == mask_token_id
            ).nonzero(as_tuple=True)[0]

            if len(mask_positions) == 0:
                logger.warning(
                    "[PepMLM] No mask positions found for length %d "
                    "(sample %d/%d); tokenizer may truncate long inputs. "
                    "Skipping.",
                    peptide_length,
                    sample_idx + 1,
                    num_peptides,
                )
                break

            # Single forward pass
            outputs = model(**inputs)
            logits = outputs.logits  # [1, seq_len, vocab_size]

            # Sample each masked position
            peptide_tokens: list[str] = []
            for pos in mask_positions:
                pos_logits = logits[0, pos, :]  # [vocab_size]

                # Restrict to amino-acid tokens
                filtered = torch.full_like(pos_logits, float("-inf"))
                for tid in aa_token_ids:
                    if tid < pos_logits.shape[0]:
                        filtered[tid] = pos_logits[tid]

                # Top-k sampling
                k = min(top_k, len(aa_token_ids))
                top_k_logits, top_k_indices = torch.topk(filtered, k)
                probs = torch.softmax(top_k_logits, dim=-1)
                chosen_local = torch.multinomial(probs, num_samples=1).item()
                chosen_token_id = top_k_indices[chosen_local].item()

                token_str = tokenizer.decode(
                    [chosen_token_id], skip_special_tokens=True
                ).strip()
                peptide_tokens.append(token_str if token_str else "G")

            peptide_seq = "".join(peptide_tokens).replace(" ", "")
            # Ensure exact length: trim or right-pad with Glycine
            peptide_seq = peptide_seq[:peptide_length]
            if len(peptide_seq) < peptide_length:
                peptide_seq = peptide_seq.ljust(peptide_length, "G")

            # Score the generated sequence
            ppl: float = float("nan")
            try:
                ppl = compute_pseudo_perplexity(
                    model, tokenizer, device, target_sequence, peptide_seq
                )
            except Exception as exc:  # noqa: BLE001
                logger.debug(
                    "[PepMLM] Perplexity scoring failed for sample %d: %s",
                    sample_idx + 1,
                    exc,
                )

            candidates.append({"sequence": peptide_seq, "perplexity": ppl})

    logger.info(
        "[PepMLM] Generated %d candidates at length %d.",
        len(candidates),
        peptide_length,
    )
    return candidates


# ---------------------------------------------------------------------------
# Pseudo-perplexity scoring
# ---------------------------------------------------------------------------

def compute_pseudo_perplexity(
    model: Any,
    tokenizer: Any,
    device: str,
    target_seq: str,
    peptide_seq: str,
) -> float:
    """Score a peptide sequence by computing its pseudo-perplexity (PPL).

    PPL = exp((1/L) * sum_i [-log P(token_i | context without token_i)])

    Each residue position in *peptide_seq* is masked one at a time; the model
    predicts the masked token from context; the negative log-likelihood of the
    true token is accumulated.  Lower PPL indicates the model assigns higher
    probability to the sequence.

    The full concatenated input (target + sep + peptide) is used so that the
    model's peptide predictions are conditioned on the target.

    Args:
        model:      Loaded PepMLM model in eval mode.
        tokenizer:  Corresponding HuggingFace tokenizer.
        device:     Torch device string.
        target_seq: Full target protein amino acid sequence.
        peptide_seq: Peptide sequence to score.

    Returns:
        Pseudo-perplexity as a float.  Returns ``inf`` if the input contains
        no scoreable positions.
    """
    torch = _import_torch()

    mask_token_id = _get_mask_token_id(tokenizer)
    sep = tokenizer.eos_token if tokenizer.eos_token else " "

    # Tokenize the full concatenation to give the model binding context
    concat_input = target_seq + sep + peptide_seq
    base_inputs = tokenizer(
        concat_input,
        return_tensors="pt",
        truncation=True,
        max_length=tokenizer.model_max_length,
    )
    full_ids = base_inputs["input_ids"].to(device)   # [1, total_len]
    attention_mask = base_inputs["attention_mask"].to(device)

    # Identify positions that correspond to the peptide portion.
    # Tokenize target alone to find where the peptide tokens begin.
    target_only = tokenizer(
        target_seq + sep,
        return_tensors="pt",
        truncation=True,
        max_length=tokenizer.model_max_length,
    )
    target_len = target_only["input_ids"].shape[1]
    total_len = full_ids.shape[1]

    # Peptide token positions: from target_len to total_len-1 (exclude final SEP)
    peptide_positions = list(range(target_len, total_len - 1))

    if not peptide_positions:
        logger.debug(
            "[PepMLM] No peptide positions available for scoring (input may be truncated)."
        )
        return float("inf")

    total_nll = 0.0
    with torch.no_grad():
        for pos in peptide_positions:
            masked_ids = full_ids.clone()
            true_token_id = int(masked_ids[0, pos].item())
            masked_ids[0, pos] = mask_token_id

            outputs = model(input_ids=masked_ids, attention_mask=attention_mask)
            logits = outputs.logits[0, pos, :]   # [vocab_size]
            log_probs = torch.log_softmax(logits, dim=-1)
            total_nll += -log_probs[true_token_id].item()

    ppl = math.exp(total_nll / len(peptide_positions))
    return ppl
