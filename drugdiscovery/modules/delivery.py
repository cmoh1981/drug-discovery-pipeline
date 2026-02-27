"""M8: Delivery system recommendation module.

Rule-based tissue-to-delivery mapping for both small molecules and peptides.
"""

from __future__ import annotations

import logging
from pathlib import Path

from drugdiscovery.types import Candidate, DeliveryRecommendation, PipelineConfig, TargetProfile
from drugdiscovery.utils.io import write_csv

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Delivery system knowledge base
# ---------------------------------------------------------------------------

SM_DELIVERY = {
    "lung": {
        "primary": "Inhalation (DPI/MDI)",
        "secondary": "IV infusion",
        "route": "inhalation",
        "rationale": "Direct pulmonary delivery maximizes local concentration and minimizes systemic exposure",
    },
    "liver": {
        "primary": "Oral (hepatic first-pass)",
        "secondary": "IV with liver-targeting nanoparticle",
        "route": "oral",
        "rationale": "Oral delivery naturally concentrates in liver via first-pass metabolism",
    },
    "brain": {
        "primary": "Intranasal (nose-to-brain)",
        "secondary": "IV with BBB shuttle peptide",
        "route": "intranasal",
        "rationale": "Intranasal bypasses BBB; IV requires BBB shuttle or nanoparticle carrier",
    },
    "blood": {
        "primary": "IV infusion",
        "secondary": "Subcutaneous",
        "route": "IV",
        "rationale": "Direct vascular delivery for blood-borne targets",
    },
    "kidney": {
        "primary": "IV (low MW for renal filtration)",
        "secondary": "Oral",
        "route": "IV",
        "rationale": "Low MW compounds naturally accumulate in kidney via glomerular filtration",
    },
    "skin": {
        "primary": "Topical (cream/gel)",
        "secondary": "Transdermal patch",
        "route": "topical",
        "rationale": "Direct dermal application for skin targets",
    },
    "eye": {
        "primary": "Intravitreal injection",
        "secondary": "Topical (eye drops)",
        "route": "intravitreal",
        "rationale": "Intravitreal for posterior segment; topical for anterior segment",
    },
    "systemic": {
        "primary": "Oral",
        "secondary": "IV infusion",
        "route": "oral",
        "rationale": "Oral preferred for systemic targets with adequate bioavailability",
    },
}

PEP_DELIVERY = {
    "lung": {
        "primary": "Inhalation (nebulized LNP)",
        "secondary": "IV with lung-targeting LNP",
        "route": "inhalation",
        "rationale": "LNP encapsulation protects peptide during nebulization; lung-targeting ligands enhance specificity",
    },
    "liver": {
        "primary": "IV with GalNAc-LNP (hepatocyte targeting)",
        "secondary": "Subcutaneous with liver-tropic LNP",
        "route": "IV",
        "rationale": "GalNAc-conjugated LNPs achieve >95% hepatocyte uptake via ASGPR",
    },
    "brain": {
        "primary": "Intranasal with CPP conjugation",
        "secondary": "IV with transferrin-receptor shuttle",
        "route": "intranasal",
        "rationale": "CPP-conjugated peptides enhance nose-to-brain transport; avoid BBB limitation",
    },
    "blood": {
        "primary": "IV infusion (PEGylated)",
        "secondary": "Subcutaneous depot",
        "route": "IV",
        "rationale": "PEGylation extends half-life for systemic peptide circulation",
    },
    "kidney": {
        "primary": "IV (LMWP carrier)",
        "secondary": "Subcutaneous",
        "route": "IV",
        "rationale": "Low-MW protamine carriers enhance renal accumulation",
    },
    "skin": {
        "primary": "Topical with CPP enhancer",
        "secondary": "Microneedle patch",
        "route": "topical",
        "rationale": "CPP enhancers improve dermal penetration of peptides",
    },
    "muscle": {
        "primary": "Intramuscular injection",
        "secondary": "Subcutaneous depot",
        "route": "intramuscular",
        "rationale": "Direct IM delivery for muscle targets",
    },
    "systemic": {
        "primary": "Subcutaneous (PEGylated)",
        "secondary": "IV infusion",
        "route": "subcutaneous",
        "rationale": "SC injection with PEGylation for sustained systemic peptide exposure",
    },
}

# Mitochondria-specific override
MITO_DELIVERY_PEP = {
    "primary": "IV with MPP conjugation (SS-31/FxR motif)",
    "secondary": "Subcutaneous with TPP-LNP",
    "route": "IV",
    "rationale": "MPP sequences (SS-31, FxR) drive mitochondrial accumulation; TPP-LNP for lipophilic carriers",
}

MITO_DELIVERY_SM = {
    "primary": "IV with TPP conjugation",
    "secondary": "Oral with mitochondrial targeting moiety",
    "route": "IV",
    "rationale": "Triphenylphosphonium (TPP) conjugation drives mitochondrial uptake via membrane potential",
}


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def recommend_delivery(
    cfg: PipelineConfig,
    candidates: list[Candidate],
    target_profile: TargetProfile,
    output_dir: Path,
) -> list[Candidate]:
    """Assign delivery system recommendations to all candidates."""
    output_dir.mkdir(parents=True, exist_ok=True)

    tissue = cfg.tissue or target_profile.target_tissue or "systemic"
    is_mito = _is_mitochondrial(target_profile)

    recommendations: list[dict] = []

    for cand in candidates:
        rec = _recommend_for_candidate(cand, tissue, is_mito)
        cand.delivery_system = rec.primary_system
        recommendations.append({
            "candidate_id": rec.candidate_id,
            "primary_system": rec.primary_system,
            "secondary_system": rec.secondary_system,
            "route": rec.route,
            "tissue": rec.tissue,
            "rationale": rec.rationale,
            "formulation_notes": rec.formulation_notes,
        })

    write_csv(output_dir / "delivery_recommendations.csv", recommendations)
    logger.info("Delivery recommendations assigned for %d candidates (tissue=%s, mito=%s)",
                len(candidates), tissue, is_mito)
    return candidates


def _recommend_for_candidate(
    cand: Candidate,
    tissue: str,
    is_mito: bool,
) -> DeliveryRecommendation:
    """Select delivery system for a single candidate."""
    tissue_lower = tissue.lower()

    if cand.modality == "peptide":
        if is_mito:
            info = MITO_DELIVERY_PEP
        else:
            info = PEP_DELIVERY.get(tissue_lower, PEP_DELIVERY["systemic"])
    else:
        if is_mito:
            info = MITO_DELIVERY_SM
        else:
            info = SM_DELIVERY.get(tissue_lower, SM_DELIVERY["systemic"])

    # Formulation notes based on candidate properties
    notes = _formulation_notes(cand, tissue_lower)

    return DeliveryRecommendation(
        candidate_id=cand.candidate_id,
        primary_system=info["primary"],
        secondary_system=info["secondary"],
        route=info["route"],
        tissue=tissue_lower,
        rationale=info["rationale"],
        formulation_notes=notes,
    )


def _is_mitochondrial(target: TargetProfile) -> bool:
    """Check if target is mitochondrial."""
    loc = target.subcellular_location.lower()
    if "mitochond" in loc:
        return True
    for term in target.go_terms:
        if "mitochond" in term.lower():
            return True
    return False


def _formulation_notes(cand: Candidate, tissue: str) -> str:
    """Generate formulation-specific notes."""
    notes = []

    if cand.modality == "peptide":
        mw = cand.molecular_weight
        if mw > 3000:
            notes.append("High MW peptide - consider LNP or nanoparticle encapsulation")
        if cand.net_charge > 3:
            notes.append("High positive charge - may enhance cell uptake but monitor hemolysis")
        if cand.gravy > 0.5:
            notes.append("Hydrophobic peptide - consider lipid-based formulation")
        if cand.modification == "pegylation":
            notes.append("PEGylated - extended half-life, may reduce immunogenicity")
        if cand.modification == "mpp_conjugation":
            notes.append("MPP-conjugated - inherent mitochondrial targeting")
        if cand.modification == "cyclization" or cand.modification == "cyclization_disulfide":
            notes.append("Cyclic - enhanced protease resistance and oral potential")
    else:
        if cand.molecular_weight > 500:
            notes.append("MW > 500 Da - reduced oral bioavailability, consider parenteral")
        if cand.net_charge < -2:
            notes.append("Highly anionic - may need salt form for formulation")

    return "; ".join(notes) if notes else "Standard formulation"
