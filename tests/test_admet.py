"""Tests for drugdiscovery.modules.admet module."""

import pytest
from drugdiscovery.modules.admet import (
    _predict_peptide_admet,
    _gravy_to_solubility,
    _predict_cpp_score,
    _predict_protease_stability,
    _predict_hemolysis_risk,
    _predict_aggregation,
    _predict_bbb,
    _compute_aggregate,
)
from drugdiscovery.types import ADMETProfile, Candidate


class TestGravyToSolubility:
    def test_hydrophilic(self):
        assert _gravy_to_solubility(-3.0) == 1.0

    def test_hydrophobic(self):
        assert _gravy_to_solubility(3.0) == 0.0

    def test_neutral(self):
        sol = _gravy_to_solubility(0.0)
        assert 0.4 < sol < 0.6


class TestCPPScore:
    def test_arginine_rich(self):
        # Tat-like CPP: RRRRRRRR
        score = _predict_cpp_score("RRRRRRRR")
        assert score > 0.5

    def test_no_positive_charge(self):
        score = _predict_cpp_score("AAAAAAA")
        assert score < 0.4

    def test_empty(self):
        assert _predict_cpp_score("") == 0.0


class TestProteaseStability:
    def test_no_sites(self):
        stability = _predict_protease_stability("AAAAAAA")
        assert stability > 0.8

    def test_many_sites(self):
        # K and R everywhere
        stability = _predict_protease_stability("KRKRKRKR")
        assert stability < 0.5


class TestHemolysisRisk:
    def test_low_risk(self):
        risk = _predict_hemolysis_risk("AAAAAAA")
        assert risk < 0.3

    def test_high_risk(self):
        # Positive + hydrophobic = hemolytic
        risk = _predict_hemolysis_risk("RRRRLLLLFFFF")
        assert risk > 0.3


class TestAggregation:
    def test_no_aggregation(self):
        score = _predict_aggregation("ARNDCE")
        assert score < 0.3

    def test_hydrophobic_stretch(self):
        score = _predict_aggregation("AAVVVVVAA")
        assert score > 0.3


class TestBBB:
    def test_small_neutral(self):
        # Small hydrophobic peptide
        score = _predict_bbb("AAVL")
        assert score > 0.3

    def test_large_charged(self):
        score = _predict_bbb("RRRRDDDDKKKKEEEEE")
        assert score < 0.3


class TestPeptideADMET:
    def test_returns_profile(self):
        cand = Candidate(
            candidate_id="TEST_001",
            modality="peptide",
            sequence="ACDEFGHIKL",
        )
        profile = _predict_peptide_admet(cand)
        assert profile.candidate_id == "TEST_001"
        assert 0.0 <= profile.aggregate_score <= 1.0
        assert isinstance(profile.flags, list)

    def test_scores_in_range(self):
        cand = Candidate(
            candidate_id="TEST_002",
            modality="peptide",
            sequence="RDFKACDE",
        )
        profile = _predict_peptide_admet(cand)
        assert 0.0 <= profile.solubility <= 1.0
        assert 0.0 <= profile.cpp_score <= 1.0
        assert 0.0 <= profile.protease_stability <= 1.0
        assert 0.0 <= profile.hemolysis_risk <= 1.0


class TestComputeAggregate:
    def test_perfect_profile(self):
        p = ADMETProfile(
            candidate_id="T1",
            solubility=1.0,
            permeability=1.0,
            hemolysis_risk=0.0,
            protease_stability=1.0,
            aggregation_propensity=0.0,
            hepatotoxicity_risk=0.0,
            herg_liability=0.0,
            oral_bioavailability=1.0,
        )
        score = _compute_aggregate(p)
        assert score > 0.9

    def test_poor_profile(self):
        p = ADMETProfile(
            candidate_id="T2",
            solubility=0.0,
            permeability=0.0,
            hemolysis_risk=1.0,
            protease_stability=0.0,
            aggregation_propensity=1.0,
            hepatotoxicity_risk=1.0,
            herg_liability=1.0,
            oral_bioavailability=0.0,
        )
        score = _compute_aggregate(p)
        assert score < 0.2
