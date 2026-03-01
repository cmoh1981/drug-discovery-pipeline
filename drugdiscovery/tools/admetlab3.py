"""ADMETlab 3.0 REST API client for small molecule ADMET prediction."""

from __future__ import annotations

import logging
from typing import Any, Optional

from drugdiscovery.utils.web import post_with_retries

logger = logging.getLogger(__name__)

ADMETLAB3_URL = "https://admetlab3.scbdd.com/api/aio-screening/predict"

# Module-level flag: once the API is confirmed unreachable (SSL error,
# connection refused, DNS failure, etc.) we skip all subsequent calls
# in the same process to avoid wasting ~12s per candidate.
_API_UNREACHABLE = False


def predict_admet_admetlab3(smiles: str) -> Optional[dict[str, Any]]:
    """Submit a SMILES string to ADMETlab 3.0 and return ADMET predictions.

    Returns dict of ADMET properties or None if API call fails.
    API docs: https://admetlab3.scbdd.com/
    """
    global _API_UNREACHABLE
    if _API_UNREACHABLE:
        return None

    try:
        resp = post_with_retries(
            ADMETLAB3_URL,
            json={"smiles": smiles},
            timeout=30,
            attempts=1,
            delay=0,
        )
        data = resp.json()

        if data.get("code") != 200 and "data" not in data:
            logger.warning("ADMETlab3 error: %s", data.get("msg", "unknown"))
            return None

        results = data.get("data", {})
        if isinstance(results, list) and results:
            results = results[0]

        return _normalize_results(results)

    except Exception as exc:
        exc_str = str(exc)
        # Cache connection-level failures so we don't retry 400+ times
        if any(kw in exc_str for kw in ("SSL", "Certificate", "ConnectionError",
                                         "Connection refused", "Name or service")):
            _API_UNREACHABLE = True
            logger.warning(
                "ADMETlab3 API unreachable (SSL/connection error) — "
                "disabling for remaining candidates: %s", exc,
            )
        else:
            logger.warning("ADMETlab3 API call failed for SMILES %.50s: %s", smiles, exc)
        return None


def predict_batch(smiles_list: list[str]) -> list[Optional[dict]]:
    """Predict ADMET for multiple SMILES. Returns list aligned with input."""
    return [predict_admet_admetlab3(s) for s in smiles_list]


def _normalize_results(raw: dict) -> dict:
    """Normalize ADMETlab3 response to standard property names."""
    return {
        "caco2_permeability": _safe_float(raw.get("Caco-2", raw.get("caco2"))),
        "hia": _safe_float(raw.get("HIA", raw.get("hia"))),
        "oral_bioavailability": _safe_float(raw.get("F(20%)", raw.get("bioavailability"))),
        "bbb_penetration": _safe_float(raw.get("BBB", raw.get("bbb"))),
        "ppb": _safe_float(raw.get("PPB", raw.get("ppb"))),
        "vd": _safe_float(raw.get("VDss", raw.get("vd"))),
        "cyp1a2_inhibitor": _safe_float(raw.get("CYP1A2-inh", raw.get("cyp1a2_inh"))),
        "cyp2c9_inhibitor": _safe_float(raw.get("CYP2C9-inh", raw.get("cyp2c9_inh"))),
        "cyp2c19_inhibitor": _safe_float(raw.get("CYP2C19-inh", raw.get("cyp2c19_inh"))),
        "cyp2d6_inhibitor": _safe_float(raw.get("CYP2D6-inh", raw.get("cyp2d6_inh"))),
        "cyp3a4_inhibitor": _safe_float(raw.get("CYP3A4-inh", raw.get("cyp3a4_inh"))),
        "t_half": _safe_float(raw.get("T1/2", raw.get("half_life"))),
        "clearance": _safe_float(raw.get("CL", raw.get("clearance"))),
        "herg_blocker": _safe_float(raw.get("hERG", raw.get("herg"))),
        "ames_toxicity": _safe_float(raw.get("AMES", raw.get("ames"))),
        "hepatotoxicity": _safe_float(raw.get("DILI", raw.get("hepatotoxicity"))),
        "skin_sensitization": _safe_float(raw.get("Skin", raw.get("skin_sens"))),
        "carcinogenicity": _safe_float(raw.get("Carcinogenicity", raw.get("carcino"))),
        "ld50": _safe_float(raw.get("LD50", raw.get("ld50"))),
        "solubility": _safe_float(raw.get("LogS", raw.get("solubility"))),
    }


def _safe_float(val, default: float = 0.0) -> float:
    if val is None:
        return default
    try:
        return float(val)
    except (TypeError, ValueError):
        return default
