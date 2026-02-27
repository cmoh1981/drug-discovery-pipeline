"""Tests for drugdiscovery.modules.scoring module."""

import pytest
from drugdiscovery.types import Candidate, TargetProfile, BindingPocket, PipelineConfig
from drugdiscovery.modules.scoring import (
    score_physicochemical,
    score_drug_likeness,
    score_selectivity,
    compute_sequence_diversity,
    compute_composite_score,
)


def _make_peptide(seq: str, cid: str = "T1") -> Candidate:
    return Candidate(
        candidate_id=cid,
        modality="peptide",
        sequence=seq,
    )


def _make_target() -> TargetProfile:
    return TargetProfile(
        gene_name="YARS2",
        uniprot_id="Q9Y2Z4",
        binding_pockets=[
            BindingPocket(
                pocket_id="P1",
                residue_numbers=[77, 81, 121],
                residue_names=["His", "Gly", "Asp"],
            ),
        ],
    )


class TestScorePhysicochemical:
    def test_peptide_mw(self):
        c = _make_peptide("ACDEFGHIKL")
        c = score_physicochemical(c)
        assert c.molecular_weight > 500
        assert c.molecular_weight < 2000

    def test_peptide_charge(self):
        c = _make_peptide("RRRR")
        c = score_physicochemical(c)
        assert c.net_charge > 2

    def test_peptide_gravy(self):
        c = _make_peptide("VVVVV")
        c = score_physicochemical(c)
        assert c.gravy > 3.0


class TestScoreDrugLikeness:
    def test_good_peptide(self):
        c = _make_peptide("ACDEFGHIKL")
        c = score_drug_likeness(c)
        assert 0.0 <= c.drug_likeness <= 1.0
        assert c.drug_likeness > 0.3

    def test_bad_peptide_too_short(self):
        c = _make_peptide("ACE")
        c = score_drug_likeness(c)
        assert c.drug_likeness < 1.0


class TestScoreSelectivity:
    def test_returns_score(self):
        c = _make_peptide("ACDEFGHIKL")
        target = _make_target()
        c = score_selectivity(c, target)
        assert 0.0 <= c.selectivity_score <= 1.0

    def test_sm_placeholder(self):
        c = Candidate(candidate_id="SM1", modality="small_molecule", smiles="c1ccccc1")
        target = _make_target()
        c = score_selectivity(c, target)
        assert c.selectivity_score == 0.5


class TestSequenceDiversity:
    def test_identical_sequences(self):
        cands = [_make_peptide("AAAA", f"T{i}") for i in range(3)]
        cands = compute_sequence_diversity(cands)
        assert all(c.sequence_diversity == 0.0 for c in cands)

    def test_different_sequences(self):
        cands = [
            _make_peptide("AAAA", "T1"),
            _make_peptide("RRRR", "T2"),
            _make_peptide("DDDD", "T3"),
        ]
        cands = compute_sequence_diversity(cands)
        assert all(c.sequence_diversity > 0.5 for c in cands)

    def test_single_candidate(self):
        cands = [_make_peptide("AAAA")]
        cands = compute_sequence_diversity(cands)
        assert cands[0].sequence_diversity == 1.0


class TestCompositeScore:
    def test_all_high(self):
        c = _make_peptide("ACDEFGHIKL")
        c.binding_score = -10.0  # good binding
        c.structure_confidence = 0.9
        c.selectivity_score = 0.9
        c.drug_likeness = 0.9
        c.sequence_diversity = 0.8
        weights = {
            "binding_energy": 0.30,
            "structure_confidence": 0.20,
            "selectivity": 0.20,
            "drug_likeness": 0.15,
            "sequence_diversity": 0.15,
        }
        # Need batch context for binding normalization
        c = compute_composite_score(c, weights, binding_min=-10.0, binding_max=0.0)
        assert 0.0 <= c.composite_score <= 1.0
        assert c.composite_score > 0.5
