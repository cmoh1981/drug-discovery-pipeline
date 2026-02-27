"""M4: Structure prediction module.

Small molecules: RDKit ETKDG conformer generation
Peptides: ESMFold / ColabFold for complex prediction
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from drugdiscovery.types import Candidate, PipelineConfig, TargetProfile

logger = logging.getLogger(__name__)


def predict_structures(
    cfg: PipelineConfig,
    candidates: list[Candidate],
    target_profile: TargetProfile,
    run_dir: Path,
) -> list[Candidate]:
    """Predict 3D structures for candidates where applicable."""
    peptides = [c for c in candidates if c.modality == "peptide"]
    small_mols = [c for c in candidates if c.modality == "small_molecule"]

    if small_mols:
        _generate_sm_conformers(small_mols, run_dir)

    if peptides:
        _predict_peptide_structures(peptides, cfg, target_profile, run_dir)

    return candidates


def _generate_sm_conformers(candidates: list[Candidate], run_dir: Path) -> None:
    """Generate 3D conformers for small molecules using RDKit ETKDG."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        logger.warning("RDKit not available; skipping SM conformer generation")
        return

    output_dir = run_dir / "de_novo" / "conformers"
    output_dir.mkdir(parents=True, exist_ok=True)

    generated = 0
    for cand in candidates:
        if not cand.smiles:
            continue
        mol = Chem.MolFromSmiles(cand.smiles)
        if mol is None:
            continue
        mol = Chem.AddHs(mol)
        result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        if result == 0:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            pdb_path = output_dir / f"{cand.candidate_id}.pdb"
            Chem.MolToPDBFile(mol, str(pdb_path))
            cand.metadata["conformer_pdb"] = str(pdb_path)
            generated += 1

    logger.info("Generated %d SM conformers out of %d", generated, len(candidates))


def _predict_peptide_structures(
    candidates: list[Candidate],
    cfg: PipelineConfig,
    target: TargetProfile,
    run_dir: Path,
) -> None:
    """Predict peptide structures using ESMFold API."""
    # Try ESMFold API (free, no auth needed)
    esmfold_ok = _try_esmfold(candidates, run_dir)

    if not esmfold_ok:
        logger.info("ESMFold unavailable; attempting ColabFold")
        _try_colabfold(candidates, cfg, target, run_dir)


def _try_esmfold(candidates: list[Candidate], run_dir: Path) -> bool:
    """Try ESMFold API for peptide structure prediction."""
    try:
        import requests
    except ImportError:
        return False

    output_dir = run_dir / "de_novo" / "structures"
    output_dir.mkdir(parents=True, exist_ok=True)

    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    predicted = 0

    # Only predict for top candidates to avoid API rate limits
    for cand in candidates[:20]:
        if not cand.sequence or len(cand.sequence) > 400:
            continue
        try:
            resp = requests.post(
                url,
                data=cand.sequence,
                headers={"Content-Type": "text/plain"},
                timeout=120,
            )
            if resp.status_code == 200:
                pdb_path = output_dir / f"{cand.candidate_id}.pdb"
                pdb_path.write_text(resp.text, encoding="utf-8")
                cand.metadata["structure_pdb"] = str(pdb_path)

                # Extract pLDDT as confidence
                plddt = _extract_mean_plddt(resp.text)
                if plddt > 0:
                    cand.structure_confidence = plddt / 100.0

                predicted += 1
            else:
                logger.debug("ESMFold failed for %s: HTTP %d", cand.candidate_id, resp.status_code)
        except Exception as exc:
            logger.debug("ESMFold error for %s: %s", cand.candidate_id, exc)
            if predicted == 0:
                return False  # API likely down

    if predicted > 0:
        logger.info("ESMFold predicted %d peptide structures", predicted)
        return True
    return False


def _try_colabfold(
    candidates: list[Candidate],
    cfg: PipelineConfig,
    target: TargetProfile,
    run_dir: Path,
) -> None:
    """Try ColabFold for complex prediction (requires local installation)."""
    try:
        from drugdiscovery.tools.colabfold import predict_complex
    except ImportError:
        logger.warning("ColabFold not available; skipping complex prediction")
        return

    output_dir = run_dir / "de_novo" / "complexes"
    output_dir.mkdir(parents=True, exist_ok=True)

    predicted = 0
    for cand in candidates[:10]:
        if not cand.sequence:
            continue
        try:
            result = predict_complex(
                target_sequence=target.sequence,
                peptide_sequence=cand.sequence,
                output_dir=str(output_dir),
                job_name=cand.candidate_id,
            )
            if result:
                cand.structure_confidence = result.get("iptm", 0.0)
                cand.metadata["complex_pdb"] = result.get("pdb_path", "")
                predicted += 1
        except Exception as exc:
            logger.debug("ColabFold failed for %s: %s", cand.candidate_id, exc)

    logger.info("ColabFold predicted %d complexes", predicted)


def _extract_mean_plddt(pdb_text: str) -> float:
    """Extract mean pLDDT from B-factor column of PDB ATOM records."""
    bfactors = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                bf = float(line[60:66].strip())
                bfactors.append(bf)
            except (ValueError, IndexError):
                pass
    if bfactors:
        return sum(bfactors) / len(bfactors)
    return 0.0
