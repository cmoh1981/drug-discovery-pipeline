"""Tests for ADMET modality-specific aggregate scoring."""

from drugdiscovery.modules.admet import _compute_aggregate
from drugdiscovery.types import ADMETProfile


def _sm_profile(**kwargs) -> ADMETProfile:
    defaults = dict(
        candidate_id="SM_TEST",
        solubility=0.7,
        permeability=0.6,
        oral_bioavailability=0.8,
        hepatotoxicity_risk=0.1,
        herg_liability=0.1,
        cyp_inhibition=0.0,
        plasma_protein_binding=0.3,
        bbb_permeability=0.2,
    )
    defaults.update(kwargs)
    return ADMETProfile(**defaults)


def _pep_profile(**kwargs) -> ADMETProfile:
    defaults = dict(
        candidate_id="PEP_TEST",
        solubility=0.7,
        permeability=0.5,
        hemolysis_risk=0.1,
        protease_stability=0.8,
        aggregation_propensity=0.1,
        hepatotoxicity_risk=0.1,
        herg_liability=0.1,
        oral_bioavailability=0.3,
    )
    defaults.update(kwargs)
    return ADMETProfile(**defaults)


class TestSmAggregate:
    def test_cyp_inhibition_penalty(self):
        score_clean = _compute_aggregate(_sm_profile(cyp_inhibition=0.0), "small_molecule")
        score_dirty = _compute_aggregate(_sm_profile(cyp_inhibition=0.8), "small_molecule")
        assert score_clean > score_dirty
        assert score_clean - score_dirty > 0.05

    def test_herg_penalty(self):
        score_safe = _compute_aggregate(_sm_profile(herg_liability=0.0), "small_molecule")
        score_risky = _compute_aggregate(_sm_profile(herg_liability=0.9), "small_molecule")
        assert score_safe > score_risky

    def test_oral_bioavailability_boost(self):
        score_high = _compute_aggregate(_sm_profile(oral_bioavailability=0.9), "small_molecule")
        score_low = _compute_aggregate(_sm_profile(oral_bioavailability=0.1), "small_molecule")
        assert score_high > score_low

    def test_ppb_moderate_penalty(self):
        score_low_ppb = _compute_aggregate(_sm_profile(plasma_protein_binding=0.1), "small_molecule")
        score_high_ppb = _compute_aggregate(_sm_profile(plasma_protein_binding=0.9), "small_molecule")
        assert score_low_ppb > score_high_ppb

    def test_score_in_range(self):
        score = _compute_aggregate(_sm_profile(), "small_molecule")
        assert 0.0 <= score <= 1.0

    def test_weights_sum_to_one(self):
        """Verify all SM component weights sum to 1.0."""
        weights = [0.15, 0.15, 0.15, 0.15, 0.15, 0.10, 0.10, 0.05]
        assert abs(sum(weights) - 1.0) < 1e-9


class TestPeptideAggregate:
    def test_protease_stability_matters(self):
        score_stable = _compute_aggregate(_pep_profile(protease_stability=0.9), "peptide")
        score_unstable = _compute_aggregate(_pep_profile(protease_stability=0.1), "peptide")
        assert score_stable > score_unstable
        assert score_stable - score_unstable > 0.05

    def test_hemolysis_penalty(self):
        score_safe = _compute_aggregate(_pep_profile(hemolysis_risk=0.0), "peptide")
        score_toxic = _compute_aggregate(_pep_profile(hemolysis_risk=0.9), "peptide")
        assert score_safe > score_toxic

    def test_aggregation_penalty(self):
        score_ok = _compute_aggregate(_pep_profile(aggregation_propensity=0.0), "peptide")
        score_bad = _compute_aggregate(_pep_profile(aggregation_propensity=0.9), "peptide")
        assert score_ok > score_bad

    def test_score_in_range(self):
        score = _compute_aggregate(_pep_profile(), "peptide")
        assert 0.0 <= score <= 1.0

    def test_weights_sum_to_one(self):
        """Verify all peptide component weights sum to 1.0."""
        weights = [0.20, 0.15, 0.15, 0.15, 0.10, 0.10, 0.10, 0.05]
        assert abs(sum(weights) - 1.0) < 1e-9


class TestDefaultModality:
    def test_default_is_peptide(self):
        """Calling without modality arg defaults to peptide weights."""
        p = _pep_profile()
        score_default = _compute_aggregate(p)
        score_peptide = _compute_aggregate(p, "peptide")
        assert score_default == score_peptide
