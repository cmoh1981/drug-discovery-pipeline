"""Tests for drugdiscovery.types module."""

from drugdiscovery.types import (
    BindingPocket,
    Candidate,
    Modality,
    ModeOfAction,
    PipelineConfig,
    TargetProfile,
    ADMETProfile,
    DeliveryRecommendation,
)


class TestModality:
    def test_values(self):
        assert Modality.SMALL_MOLECULE == "small_molecule"
        assert Modality.PEPTIDE == "peptide"


class TestModeOfAction:
    def test_values(self):
        assert ModeOfAction.AGONIST == "agonist"
        assert ModeOfAction.ANTAGONIST == "antagonist"


class TestCandidate:
    def test_defaults(self):
        c = Candidate(candidate_id="TEST_001")
        assert c.candidate_id == "TEST_001"
        assert c.composite_score == 0.0
        assert c.rank == 0

    def test_to_dict(self):
        c = Candidate(
            candidate_id="TEST_001",
            modality="peptide",
            sequence="ACDEFG",
            composite_score=0.8765,
        )
        d = c.to_dict()
        assert d["candidate_id"] == "TEST_001"
        assert d["composite_score"] == 0.8765
        assert isinstance(d, dict)

    def test_to_dict_rounding(self):
        c = Candidate(
            candidate_id="T1",
            molecular_weight=1234.56789,
            binding_score=0.123456789,
        )
        d = c.to_dict()
        assert d["molecular_weight"] == 1234.57
        assert d["binding_score"] == 0.1235


class TestTargetProfile:
    def test_defaults(self):
        tp = TargetProfile(gene_name="YARS2", uniprot_id="Q9Y2Z4")
        assert tp.gene_name == "YARS2"
        assert tp.binding_pockets == []
        assert tp.anti_targets == []

    def test_binding_pockets(self):
        pocket = BindingPocket(
            pocket_id="P1",
            residue_numbers=[77, 81, 121],
            description="Active site",
        )
        tp = TargetProfile(
            gene_name="YARS2",
            uniprot_id="Q9Y2Z4",
            binding_pockets=[pocket],
        )
        assert len(tp.binding_pockets) == 1
        assert tp.binding_pockets[0].residue_numbers == [77, 81, 121]


class TestPipelineConfig:
    def test_defaults(self):
        cfg = PipelineConfig()
        assert cfg.top_n == 20
        assert cfg.device == "auto"
        assert sum(cfg.scoring_weights.values()) == 1.0

    def test_scoring_weights_sum(self):
        cfg = PipelineConfig()
        total = sum(cfg.scoring_weights.values())
        assert abs(total - 1.0) < 0.001


class TestADMETProfile:
    def test_defaults(self):
        p = ADMETProfile(candidate_id="T1")
        assert p.flag_count == 0
        assert p.flags == []


class TestDeliveryRecommendation:
    def test_creation(self):
        r = DeliveryRecommendation(
            candidate_id="T1",
            primary_system="IV infusion",
            route="IV",
        )
        assert r.primary_system == "IV infusion"
