"""M2: Library screening module.

Queries multiple databases, merges results, and deduplicates candidates.
"""

from __future__ import annotations

import logging
from pathlib import Path

from drugdiscovery.types import Candidate, Modality, PipelineConfig, TargetProfile
from drugdiscovery.utils.io import write_csv

logger = logging.getLogger(__name__)


def screen_libraries(
    cfg: PipelineConfig,
    target: TargetProfile,
    output_dir: Path,
) -> list[Candidate]:
    """Query all relevant databases and return merged candidates."""
    output_dir.mkdir(parents=True, exist_ok=True)
    all_hits: list[Candidate] = []
    limit = cfg.library_settings.get("max_hits_per_db", 500)

    if cfg.modality == Modality.SMALL_MOLECULE:
        all_hits.extend(_screen_sm_databases(target, limit, output_dir))
    else:
        all_hits.extend(_screen_peptide_databases(target, limit, output_dir))

    # Deduplicate
    before = len(all_hits)
    all_hits = _deduplicate(all_hits)
    logger.info("Deduplicated: %d -> %d candidates", before, len(all_hits))

    # Save merged
    if all_hits:
        rows = [c.to_dict() for c in all_hits]
        write_csv(output_dir / "merged_hits.csv", rows)

    return all_hits


def _screen_sm_databases(
    target: TargetProfile,
    limit: int,
    output_dir: Path,
) -> list[Candidate]:
    """Screen small molecule databases."""
    hits: list[Candidate] = []

    # ChEMBL
    try:
        from drugdiscovery.databases.chembl import search_by_target

        results = search_by_target(target.gene_name, limit=limit)
        for r in results:
            c = Candidate(
                candidate_id=r.get("id", f"CHEMBL_{len(hits)+1}"),
                candidate_type="library_hit",
                modality="small_molecule",
                source="chembl",
                smiles=r.get("smiles", ""),
                molecular_weight=_safe_float(r.get("mw", 0)),
                moa_predicted=r.get("moa", "unknown"),
            )
            if c.smiles:
                hits.append(c)
        _save_db_hits(output_dir / "chembl_hits.csv", results)
        logger.info("ChEMBL: %d hits for %s", len(results), target.gene_name)
    except Exception as exc:
        logger.warning("ChEMBL search failed: %s", exc)

    # DrugBank
    try:
        from drugdiscovery.databases.drugbank import search_drugbank

        results = search_drugbank(target.gene_name, limit=limit)
        for r in results:
            c = Candidate(
                candidate_id=r.get("id", f"DB_{len(hits)+1}"),
                candidate_type="library_hit",
                modality="small_molecule",
                source="drugbank",
                smiles=r.get("smiles", ""),
                moa_predicted=r.get("moa", "unknown"),
            )
            if c.smiles:
                hits.append(c)
        _save_db_hits(output_dir / "drugbank_hits.csv", results)
        logger.info("DrugBank: %d hits", len(results))
    except Exception as exc:
        logger.warning("DrugBank search failed: %s", exc)

    # COCONUT
    try:
        from drugdiscovery.databases.coconut import search_coconut

        results = search_coconut(target.gene_name, limit=limit)
        for r in results:
            c = Candidate(
                candidate_id=r.get("id", f"COCONUT_{len(hits)+1}"),
                candidate_type="library_hit",
                modality="small_molecule",
                source="coconut",
                smiles=r.get("smiles", ""),
                molecular_weight=_safe_float(r.get("molecular_weight", 0)),
            )
            if c.smiles:
                hits.append(c)
        _save_db_hits(output_dir / "coconut_hits.csv", results)
        logger.info("COCONUT: %d hits", len(results))
    except Exception as exc:
        logger.warning("COCONUT search failed: %s", exc)

    # ZINC22
    try:
        from drugdiscovery.databases.zinc22 import search_zinc_by_similarity

        # Use first ChEMBL hit SMILES as query if available
        query_smiles = next((c.smiles for c in hits if c.smiles), "")
        if query_smiles:
            results = search_zinc_by_similarity(query_smiles, limit=limit)
            for r in results:
                c = Candidate(
                    candidate_id=r.get("id", f"ZINC_{len(hits)+1}"),
                    candidate_type="library_hit",
                    modality="small_molecule",
                    source="zinc22",
                    smiles=r.get("smiles", ""),
                )
                if c.smiles:
                    hits.append(c)
            _save_db_hits(output_dir / "zinc22_hits.csv", results)
            logger.info("ZINC22: %d hits", len(results))
    except Exception as exc:
        logger.warning("ZINC22 search failed: %s", exc)

    return hits


def _screen_peptide_databases(
    target: TargetProfile,
    limit: int,
    output_dir: Path,
) -> list[Candidate]:
    """Screen peptide databases."""
    hits: list[Candidate] = []

    # ChEMBL peptides
    try:
        from drugdiscovery.databases.chembl import search_by_target

        results = search_by_target(target.gene_name, limit=limit)
        for r in results:
            smiles = r.get("smiles", "")
            # Filter for peptide-like molecules (MW > 500, has amino acid substructure)
            mw = _safe_float(r.get("mw", 0))
            if mw > 500:
                c = Candidate(
                    candidate_id=r.get("id", f"CHEMBL_PEP_{len(hits)+1}"),
                    candidate_type="library_hit",
                    modality="peptide",
                    source="chembl",
                    smiles=smiles,
                    molecular_weight=mw,
                )
                hits.append(c)
        logger.info("ChEMBL peptides: %d hits", len([h for h in hits if h.source == "chembl"]))
    except Exception as exc:
        logger.warning("ChEMBL peptide search failed: %s", exc)

    # BIOPEP
    try:
        from drugdiscovery.databases.peptide_dbs import search_biopep

        results = search_biopep(activity="inhibitor", limit=limit)
        for r in results:
            c = Candidate(
                candidate_id=r.get("id", f"BIOPEP_{len(hits)+1}"),
                candidate_type="library_hit",
                modality="peptide",
                source="biopep",
                sequence=r.get("sequence", ""),
            )
            if c.sequence:
                hits.append(c)
        _save_db_hits(output_dir / "biopep_hits.csv", results)
        logger.info("BIOPEP: %d hits", len(results))
    except Exception as exc:
        logger.warning("BIOPEP search failed: %s", exc)

    # APD3
    try:
        from drugdiscovery.databases.peptide_dbs import search_apd3

        results = search_apd3(limit=limit)
        for r in results:
            c = Candidate(
                candidate_id=r.get("id", f"APD3_{len(hits)+1}"),
                candidate_type="library_hit",
                modality="peptide",
                source="apd3",
                sequence=r.get("sequence", ""),
            )
            if c.sequence:
                hits.append(c)
        _save_db_hits(output_dir / "apd3_hits.csv", results)
        logger.info("APD3: %d hits", len(results))
    except Exception as exc:
        logger.warning("APD3 search failed: %s", exc)

    # DBAASP
    try:
        from drugdiscovery.databases.peptide_dbs import search_dbaasp

        results = search_dbaasp(target=target.gene_name, limit=limit)
        for r in results:
            c = Candidate(
                candidate_id=r.get("id", f"DBAASP_{len(hits)+1}"),
                candidate_type="library_hit",
                modality="peptide",
                source="dbaasp",
                sequence=r.get("sequence", ""),
            )
            if c.sequence:
                hits.append(c)
        _save_db_hits(output_dir / "dbaasp_hits.csv", results)
        logger.info("DBAASP: %d hits", len(results))
    except Exception as exc:
        logger.warning("DBAASP search failed: %s", exc)

    return hits


def _deduplicate(candidates: list[Candidate]) -> list[Candidate]:
    """Remove duplicate candidates by SMILES or sequence."""
    seen: set[str] = set()
    unique: list[Candidate] = []

    for c in candidates:
        key = c.smiles.strip() if c.smiles else c.sequence.strip()
        if not key or key in seen:
            continue
        seen.add(key)
        unique.append(c)

    return unique


def _save_db_hits(path: Path, results: list[dict]) -> None:
    """Save raw database hits to CSV."""
    if results:
        write_csv(path, results)


def _safe_float(val, default: float = 0.0) -> float:
    try:
        return float(val)
    except (TypeError, ValueError):
        return default
