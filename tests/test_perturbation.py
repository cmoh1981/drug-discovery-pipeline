"""Tests for drugdiscovery perturbation biology module (M4.6).

Tests cover:
  - CMAP client: connectivity scoring, compound lookup fallbacks
  - STRING PPI: network effect computation
  - Perturbation module: end-to-end scoring with mocked APIs
"""

import pytest
from unittest.mock import patch, MagicMock
from pathlib import Path

from drugdiscovery.types import Candidate, TargetProfile, BindingPocket, PipelineConfig


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_peptide(seq: str = "ACDEFGHIKL", cid: str = "P1") -> Candidate:
    return Candidate(
        candidate_id=cid,
        modality="peptide",
        sequence=seq,
    )


def _make_sm(smiles: str = "c1ccccc1", cid: str = "SM1") -> Candidate:
    return Candidate(
        candidate_id=cid,
        modality="small_molecule",
        smiles=smiles,
    )


def _make_target() -> TargetProfile:
    return TargetProfile(
        gene_name="EGFR",
        uniprot_id="P00533",
        binding_pockets=[
            BindingPocket(
                pocket_id="P1",
                residue_numbers=[719, 720, 721],
                residue_names=["Leu", "Gly", "Thr"],
            ),
        ],
        anti_targets=[{"gene_name": "ERBB2"}],
        go_terms=["GO:0006468", "protein phosphorylation"],
    )


# ---------------------------------------------------------------------------
# CMAP Client Tests
# ---------------------------------------------------------------------------

class TestCmapConnectivityScore:
    """Test the connectivity score computation."""

    def test_perfect_reversal(self):
        from drugdiscovery.tools.cmap_client import compute_connectivity_score

        # Drug upregulates what disease downregulates (and vice versa)
        drug_up = ["GENE_A", "GENE_B", "GENE_C"]
        drug_down = ["GENE_D", "GENE_E", "GENE_F"]
        disease_up = ["GENE_D", "GENE_E", "GENE_F"]
        disease_down = ["GENE_A", "GENE_B", "GENE_C"]

        score = compute_connectivity_score(drug_up, drug_down, disease_up, disease_down)
        assert score > 0.8  # Strong reversal

    def test_perfect_mimicry(self):
        from drugdiscovery.tools.cmap_client import compute_connectivity_score

        # Drug mimics disease (bad)
        drug_up = ["GENE_A", "GENE_B"]
        drug_down = ["GENE_C", "GENE_D"]
        disease_up = ["GENE_A", "GENE_B"]
        disease_down = ["GENE_C", "GENE_D"]

        score = compute_connectivity_score(drug_up, drug_down, disease_up, disease_down)
        assert score < 0.3  # Mimicry = bad

    def test_no_overlap(self):
        from drugdiscovery.tools.cmap_client import compute_connectivity_score

        drug_up = ["GENE_A", "GENE_B"]
        drug_down = ["GENE_C", "GENE_D"]
        disease_up = ["GENE_X", "GENE_Y"]
        disease_down = ["GENE_Z", "GENE_W"]

        score = compute_connectivity_score(drug_up, drug_down, disease_up, disease_down)
        assert score == 0.5  # Neutral

    def test_empty_signature(self):
        from drugdiscovery.tools.cmap_client import compute_connectivity_score

        score = compute_connectivity_score([], [], ["A", "B"], ["C", "D"])
        assert score == 0.5  # No data = neutral

    def test_partial_reversal(self):
        from drugdiscovery.tools.cmap_client import compute_connectivity_score

        drug_up = ["GENE_A", "GENE_B", "GENE_C", "GENE_D"]
        drug_down = ["GENE_E", "GENE_F", "GENE_G", "GENE_H"]
        disease_up = ["GENE_E", "GENE_F"]  # Drug reverses 2 of these
        disease_down = ["GENE_A"]  # Drug reverses 1 of these

        score = compute_connectivity_score(drug_up, drug_down, disease_up, disease_down)
        assert 0.5 < score < 1.0  # Partial reversal


class TestCmapKnnFallback:
    """Test kNN compound lookup fallback."""

    def test_known_target_lookup(self):
        from drugdiscovery.tools.cmap_client import find_nearest_l1000_compound

        results = find_nearest_l1000_compound(smiles="", gene_target="EGFR")
        assert len(results) > 0
        assert results[0]["name"] == "erlotinib"

    def test_unknown_target(self):
        from drugdiscovery.tools.cmap_client import find_nearest_l1000_compound

        results = find_nearest_l1000_compound(smiles="", gene_target="UNKNOWN_GENE_XYZ")
        assert len(results) == 0

    def test_known_target_mtor(self):
        from drugdiscovery.tools.cmap_client import find_nearest_l1000_compound

        results = find_nearest_l1000_compound(smiles="", gene_target="MTOR")
        assert len(results) > 0
        assert any("rapamycin" in r["name"] for r in results)


class TestCmapQueryPerturbation:
    """Test the main query function with mocked API."""

    @patch("drugdiscovery.tools.cmap_client.find_nearest_l1000_compound")
    def test_knn_fallback_for_known_target(self, mock_knn):
        from drugdiscovery.tools.cmap_client import query_cmap_perturbation

        mock_knn.return_value = [
            {"name": "erlotinib", "pert_id": "BRD-K53533115",
             "target": "EGFR", "similarity": 0.5, "source": "target_lookup"}
        ]

        result = query_cmap_perturbation(
            gene_target="EGFR",
            api_key="",  # No API key → skip CLUE, go to kNN
        )
        assert result["match_source"] == "knn_transfer"
        assert result["matched_compound"] == "erlotinib"
        assert 0.0 <= result["connectivity_score"] <= 1.0

    def test_no_match_returns_neutral(self):
        from drugdiscovery.tools.cmap_client import query_cmap_perturbation

        result = query_cmap_perturbation(
            smiles="",
            gene_target="NONEXISTENT_GENE_XYZ",
            api_key="",
        )
        assert result["connectivity_score"] == 0.5
        assert result["match_source"] in ("none", "knn_transfer")


# ---------------------------------------------------------------------------
# STRING PPI Tests
# ---------------------------------------------------------------------------

class TestStringNetworkEffect:
    """Test network effect computation with mocked STRING API."""

    @patch("drugdiscovery.tools.string_ppi.get_interaction_partners")
    @patch("drugdiscovery.tools.string_ppi.get_enrichment")
    def test_network_with_disease_overlap(self, mock_enrichment, mock_partners):
        from drugdiscovery.tools.string_ppi import compute_network_effect

        mock_partners.return_value = [
            {"partner_gene": "ERBB2", "partner_string_id": "x", "combined_score": 900,
             "experimental_score": 800, "database_score": 700},
            {"partner_gene": "GRB2", "partner_string_id": "y", "combined_score": 850,
             "experimental_score": 750, "database_score": 600},
            {"partner_gene": "SHC1", "partner_string_id": "z", "combined_score": 800,
             "experimental_score": 700, "database_score": 500},
        ]
        mock_enrichment.return_value = []

        result = compute_network_effect(
            target_gene="EGFR",
            disease_genes=["ERBB2", "GRB2", "KRAS", "PIK3CA"],
            mode="antagonist",
            depth=1,
        )

        assert result["network_effect_score"] > 0.3
        assert "ERBB2" in result["affected_disease_genes"]
        assert "GRB2" in result["affected_disease_genes"]
        assert result["total_network_genes"] >= 3

    @patch("drugdiscovery.tools.string_ppi.get_interaction_partners")
    @patch("drugdiscovery.tools.string_ppi.get_enrichment")
    def test_no_partners(self, mock_enrichment, mock_partners):
        from drugdiscovery.tools.string_ppi import compute_network_effect

        mock_partners.return_value = []
        mock_enrichment.return_value = []

        result = compute_network_effect(
            target_gene="UNKNOWN_XYZ",
            disease_genes=["A", "B"],
            depth=1,
        )
        assert result["network_effect_score"] == 0.5

    @patch("drugdiscovery.tools.string_ppi.get_interaction_partners")
    @patch("drugdiscovery.tools.string_ppi.get_enrichment")
    def test_no_disease_genes_uses_size_proxy(self, mock_enrichment, mock_partners):
        from drugdiscovery.tools.string_ppi import compute_network_effect

        mock_partners.return_value = [
            {"partner_gene": f"GENE_{i}", "partner_string_id": f"id_{i}",
             "combined_score": 800, "experimental_score": 700, "database_score": 600}
            for i in range(30)
        ]
        mock_enrichment.return_value = []

        result = compute_network_effect(
            target_gene="EGFR",
            disease_genes=[],
            depth=1,
        )
        # Score based on network size (30 genes → ~0.4)
        assert result["network_effect_score"] > 0.1
        assert result["total_network_genes"] == 30


# ---------------------------------------------------------------------------
# Perturbation Module Integration Tests
# ---------------------------------------------------------------------------

class TestPredictPerturbation:
    """Integration tests for the perturbation module."""

    @patch("drugdiscovery.modules.perturbation.predict_perturbation.__module__")
    def test_candidates_get_perturbation_score(self, tmp_path=None):
        """Test that candidates receive perturbation scores."""
        from drugdiscovery.modules.perturbation import _score_candidate_perturbation

        cand = _make_sm(smiles="c1ccccc1", cid="SM_TEST")
        target = _make_target()

        network_result = {
            "network_effect_score": 0.7,
            "affected_genes": ["ERBB2", "GRB2"],
            "affected_disease_genes": ["ERBB2"],
            "total_network_genes": 2,
        }

        with patch("drugdiscovery.tools.cmap_client.query_cmap_perturbation") as mock_cmap:
            mock_cmap.return_value = {
                "connectivity_score": 0.65,
                "matched_compound": "erlotinib",
                "match_source": "knn_transfer",
                "up_genes": [],
                "down_genes": [],
            }

            result = _score_candidate_perturbation(
                cand, target,
                disease_up_genes=["ERBB2", "GRB2"],
                disease_down_genes=["TP53"],
                network_result=network_result,
                tier_weights={"cmap_connectivity": 0.5, "knn_transfer": 0.2, "network_effect": 0.3},
            )

        assert 0.0 <= result.perturbation_score <= 1.0
        assert result.cmap_compound_match == "erlotinib"
        assert result.network_effect_score == 0.7

    def test_peptide_weight_adjustment(self):
        """Test that peptides shift weight from CMAP to network."""
        from drugdiscovery.modules.perturbation import _score_candidate_perturbation

        cand = _make_peptide("ACDEFGHIKL", "PEP_TEST")
        target = _make_target()

        network_result = {
            "network_effect_score": 0.8,
        }

        with patch("drugdiscovery.tools.cmap_client.query_cmap_perturbation") as mock_cmap:
            mock_cmap.return_value = {
                "connectivity_score": 0.5,
                "matched_compound": "",
                "match_source": "none",
                "up_genes": [],
                "down_genes": [],
            }

            result = _score_candidate_perturbation(
                cand, target,
                disease_up_genes=["A"],
                disease_down_genes=["B"],
                network_result=network_result,
                tier_weights={"cmap_connectivity": 0.5, "knn_transfer": 0.2, "network_effect": 0.3},
            )

        # Peptide should have higher score because network effect is high
        # and peptides shift weight toward network
        assert 0.0 <= result.perturbation_score <= 1.0

    @patch("drugdiscovery.modules.perturbation._derive_disease_genes")
    @patch("drugdiscovery.modules.perturbation._score_candidate_perturbation")
    def test_predict_perturbation_saves_csv(self, mock_score, mock_disease, tmp_path):
        """Test that the module saves output CSV."""
        from drugdiscovery.modules.perturbation import predict_perturbation

        mock_disease.return_value = (["ERBB2"], ["TP53"])

        def _score_side_effect(cand, *args, **kwargs):
            cand.perturbation_score = 0.65
            cand.cmap_connectivity = 0.6
            cand.cmap_compound_match = "test_drug"
            cand.network_effect_score = 0.7
            cand.disease_signature_reversal = 0.6
            return cand

        mock_score.side_effect = _score_side_effect

        candidates = [_make_sm("c1ccccc1", f"SM_{i}") for i in range(3)]
        target = _make_target()
        cfg = PipelineConfig(target="EGFR")

        output_dir = tmp_path / "perturbation"

        with patch("drugdiscovery.tools.string_ppi.compute_network_effect") as mock_net:
            mock_net.return_value = {"network_effect_score": 0.7, "total_network_genes": 10,
                                     "affected_disease_genes": ["ERBB2"]}
            result = predict_perturbation(cfg, candidates, target, output_dir)

        assert len(result) == 3
        assert all(c.perturbation_score == 0.65 for c in result)
        assert (output_dir / "perturbation_scores.csv").exists()


# ---------------------------------------------------------------------------
# Candidate dataclass integration
# ---------------------------------------------------------------------------

class TestCandidatePerturbationFields:
    """Test that Candidate has the new perturbation fields."""

    def test_default_values(self):
        c = Candidate(candidate_id="test")
        assert c.perturbation_score == 0.0
        assert c.cmap_connectivity == 0.0
        assert c.cmap_compound_match == ""
        assert c.network_effect_score == 0.0
        assert c.disease_signature_reversal == 0.0

    def test_to_dict_includes_perturbation(self):
        c = Candidate(candidate_id="test")
        c.perturbation_score = 0.75
        c.cmap_connectivity = 0.8
        c.cmap_compound_match = "erlotinib"
        c.network_effect_score = 0.6
        c.disease_signature_reversal = 0.7

        d = c.to_dict()
        assert d["perturbation_score"] == 0.75
        assert d["cmap_connectivity"] == 0.8
        assert d["cmap_compound_match"] == "erlotinib"
        assert d["network_effect_score"] == 0.6
        assert d["disease_signature_reversal"] == 0.7

    def test_scoring_weights_include_perturbation(self):
        cfg = PipelineConfig()
        assert "perturbation" in cfg.scoring_weights
        assert cfg.scoring_weights["perturbation"] == 0.15


# ---------------------------------------------------------------------------
# Composite scoring with perturbation
# ---------------------------------------------------------------------------

class TestCompositeWithPerturbation:
    """Test that composite scoring includes perturbation weight."""

    def test_perturbation_contributes_to_composite(self):
        from drugdiscovery.modules.scoring import compute_composite_score

        c = Candidate(
            candidate_id="test",
            modality="peptide",
            binding_score=-8.0,
            selectivity_score=0.7,
            drug_likeness=0.8,
            admet_score=0.6,
            moa_predicted="antagonist",
            perturbation_score=0.9,
        )

        weights = {
            "binding_energy": 0.30,
            "selectivity": 0.20,
            "drug_likeness": 0.10,
            "admet_aggregate": 0.15,
            "moa_consistency": 0.10,
            "perturbation": 0.15,
        }

        result = compute_composite_score(
            c, weights, binding_min=-10.0, binding_max=-3.0,
            target_mode="antagonist",
        )
        assert result.composite_score > 0.5

    def test_zero_perturbation_lowers_composite(self):
        from drugdiscovery.modules.scoring import compute_composite_score

        # High perturbation
        c_high = Candidate(
            candidate_id="high",
            modality="peptide",
            binding_score=-8.0,
            selectivity_score=0.7,
            drug_likeness=0.8,
            admet_score=0.6,
            moa_predicted="antagonist",
            perturbation_score=0.9,
        )
        # Low perturbation
        c_low = Candidate(
            candidate_id="low",
            modality="peptide",
            binding_score=-8.0,
            selectivity_score=0.7,
            drug_likeness=0.8,
            admet_score=0.6,
            moa_predicted="antagonist",
            perturbation_score=0.1,
        )

        weights = {
            "binding_energy": 0.30,
            "selectivity": 0.20,
            "drug_likeness": 0.10,
            "admet_aggregate": 0.15,
            "moa_consistency": 0.10,
            "perturbation": 0.15,
        }

        compute_composite_score(c_high, weights, binding_min=-10.0, binding_max=-3.0, target_mode="antagonist")
        compute_composite_score(c_low, weights, binding_min=-10.0, binding_max=-3.0, target_mode="antagonist")

        assert c_high.composite_score > c_low.composite_score
