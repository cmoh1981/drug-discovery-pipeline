"""Tests for drugdiscovery.modules.de_novo_sm module."""

import tempfile
from pathlib import Path
from unittest.mock import patch

from drugdiscovery.modules.de_novo_sm import (
    _brics_recombination,
    _rdkit_fragment_enumeration,
    _tanimoto_diversity_filter,
    _try_reinvent4,
    design_small_molecules,
)
from drugdiscovery.types import Candidate, PipelineConfig, TargetProfile


def _make_cfg(**overrides) -> PipelineConfig:
    cfg = PipelineConfig()
    for k, v in overrides.items():
        setattr(cfg, k, v)
    return cfg


def _make_target() -> TargetProfile:
    return TargetProfile(gene_name="EGFR", uniprot_id="P00533")


# ---------------------------------------------------------------------------
# _try_reinvent4  (always returns [] when reinvent not installed)
# ---------------------------------------------------------------------------


class TestTryReinvent4:
    def test_returns_empty_when_not_installed(self):
        cfg = _make_cfg()
        target = _make_target()
        with tempfile.TemporaryDirectory() as d:
            result = _try_reinvent4(cfg, target, Path(d))
        assert result == []


# ---------------------------------------------------------------------------
# _brics_recombination
# ---------------------------------------------------------------------------


class TestBricsRecombination:
    def test_generates_candidates(self):
        cfg = _make_cfg()
        target = _make_target()
        with tempfile.TemporaryDirectory() as d:
            candidates = _brics_recombination(cfg, target, Path(d))
        assert len(candidates) > 0
        assert all(c.source == "brics_recombination" for c in candidates)
        assert all(c.modality == "small_molecule" for c in candidates)
        assert all(c.candidate_type == "de_novo" for c in candidates)

    def test_all_have_valid_smiles(self):
        from rdkit import Chem

        cfg = _make_cfg()
        target = _make_target()
        with tempfile.TemporaryDirectory() as d:
            candidates = _brics_recombination(cfg, target, Path(d))
        for c in candidates:
            assert c.smiles, f"{c.candidate_id} has no SMILES"
            mol = Chem.MolFromSmiles(c.smiles)
            assert mol is not None, f"{c.candidate_id} has invalid SMILES: {c.smiles}"

    def test_molecular_weight_in_drug_range(self):
        cfg = _make_cfg()
        target = _make_target()
        with tempfile.TemporaryDirectory() as d:
            candidates = _brics_recombination(cfg, target, Path(d))
        for c in candidates:
            assert 150 < c.molecular_weight < 600, (
                f"{c.candidate_id} MW={c.molecular_weight} out of range"
            )

    def test_all_unique_smiles(self):
        cfg = _make_cfg()
        target = _make_target()
        with tempfile.TemporaryDirectory() as d:
            candidates = _brics_recombination(cfg, target, Path(d))
        smiles = [c.smiles for c in candidates]
        assert len(smiles) == len(set(smiles)), "Duplicate SMILES found"

    def test_seed_smiles_contribute_fragments(self):
        """Seed SMILES should produce additional fragments and candidates."""
        cfg = _make_cfg()
        target = _make_target()
        seeds = ["c1ccc(NC(=O)CCl)cc1", "O=C(O)c1cccc(O)c1O"]
        with tempfile.TemporaryDirectory() as d:
            candidates = _brics_recombination(cfg, target, Path(d), seed_smiles=seeds)
        assert len(candidates) > 0

    def test_candidate_ids_have_brics_prefix(self):
        cfg = _make_cfg()
        target = _make_target()
        with tempfile.TemporaryDirectory() as d:
            candidates = _brics_recombination(cfg, target, Path(d))
        for c in candidates:
            assert c.candidate_id.startswith("BRICS_")


# ---------------------------------------------------------------------------
# _tanimoto_diversity_filter
# ---------------------------------------------------------------------------


class TestTanimotoDiversityFilter:
    def test_removes_similar_molecules(self):
        """Identical SMILES should be filtered to just one."""
        candidates = [
            Candidate(candidate_id="A", smiles="c1ccccc1"),
            Candidate(candidate_id="B", smiles="c1ccccc1"),  # same as A
            Candidate(candidate_id="C", smiles="c1ccncc1"),  # different
        ]
        result = _tanimoto_diversity_filter(candidates, threshold=0.85)
        # A and B are identical (Tanimoto=1.0), so B should be removed
        assert len(result) <= 2
        assert result[0].candidate_id == "A"

    def test_keeps_diverse_molecules(self):
        candidates = [
            Candidate(candidate_id="A", smiles="c1ccccc1"),           # benzene
            Candidate(candidate_id="B", smiles="C1CCCCC1"),           # cyclohexane
            Candidate(candidate_id="C", smiles="c1ccncc1"),           # pyridine
            Candidate(candidate_id="D", smiles="c1ccc2[nH]ccc2c1"),  # indole
        ]
        result = _tanimoto_diversity_filter(candidates, threshold=0.85)
        assert len(result) >= 3  # these are all fairly different

    def test_empty_input(self):
        result = _tanimoto_diversity_filter([], threshold=0.85)
        assert result == []

    def test_single_candidate(self):
        candidates = [Candidate(candidate_id="A", smiles="c1ccccc1")]
        result = _tanimoto_diversity_filter(candidates, threshold=0.85)
        assert len(result) == 1


# ---------------------------------------------------------------------------
# _rdkit_fragment_enumeration
# ---------------------------------------------------------------------------


class TestRdkitFragmentEnumeration:
    def test_generates_candidates(self):
        cfg = _make_cfg()
        target = _make_target()
        with tempfile.TemporaryDirectory() as d:
            candidates = _rdkit_fragment_enumeration(cfg, target, Path(d))
        assert len(candidates) > 0
        assert all(c.source == "rdkit_enum" for c in candidates)

    def test_all_unique(self):
        cfg = _make_cfg()
        target = _make_target()
        with tempfile.TemporaryDirectory() as d:
            candidates = _rdkit_fragment_enumeration(cfg, target, Path(d))
        smiles = [c.smiles for c in candidates]
        assert len(smiles) == len(set(smiles))


# ---------------------------------------------------------------------------
# design_small_molecules (integration)
# ---------------------------------------------------------------------------


class TestDesignSmallMolecules:
    def test_produces_candidates_and_csv(self):
        cfg = _make_cfg()
        target = _make_target()
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "de_novo"
            candidates = design_small_molecules(cfg, target, out)
            assert len(candidates) > 0
            assert (out / "candidates.csv").exists()

    def test_falls_through_to_rdkit_when_brics_fails(self):
        """When BRICS returns nothing, should fall through to RDKit enum."""
        cfg = _make_cfg()
        target = _make_target()
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "de_novo"
            with patch(
                "drugdiscovery.modules.de_novo_sm._brics_recombination",
                return_value=[],
            ):
                candidates = design_small_molecules(cfg, target, out)
            assert len(candidates) > 0
            assert all(c.source == "rdkit_enum" for c in candidates)
