"""Tests for drugdiscovery.utils.chemistry module."""

import pytest
from drugdiscovery.utils.chemistry import (
    compute_molecular_weight,
    compute_net_charge,
    compute_gravy,
    compute_isoelectric_point,
    compute_instability_index,
    compute_aromaticity,
    count_protease_sites,
    check_drug_likeness,
    compute_peptide_properties,
    _clean_sequence,
    safe_float,
)


class TestCleanSequence:
    def test_uppercase(self):
        assert _clean_sequence("acdefg") == "ACDEFG"

    def test_strip_whitespace(self):
        assert _clean_sequence("  ACDE  ") == "ACDE"

    def test_remove_non_standard(self):
        assert _clean_sequence("ACXBZOU") == "AC"

    def test_empty(self):
        assert _clean_sequence("") == ""


class TestMolecularWeight:
    def test_single_residue(self):
        mw = compute_molecular_weight("A")
        assert mw == pytest.approx(89.09, abs=0.1)

    def test_dipeptide(self):
        # Two alanines minus one water
        mw = compute_molecular_weight("AA")
        expected = 89.09 * 2 - 18.015
        assert mw == pytest.approx(expected, abs=0.1)

    def test_known_peptide(self):
        # RDFK (SS-31 motif)
        mw = compute_molecular_weight("RDFK")
        assert 500 < mw < 600

    def test_empty(self):
        assert compute_molecular_weight("") == 0.0


class TestNetCharge:
    def test_neutral(self):
        # Alanine-only peptide should be roughly neutral
        charge = compute_net_charge("AAAA")
        assert -1 < charge < 1

    def test_positive(self):
        # Arginine-rich should be positive
        charge = compute_net_charge("RRRR")
        assert charge > 2

    def test_negative(self):
        # Aspartate-rich should be negative
        charge = compute_net_charge("DDDD")
        assert charge < -2


class TestGRAVY:
    def test_hydrophobic(self):
        gravy = compute_gravy("VVVVV")
        assert gravy > 3.0

    def test_hydrophilic(self):
        gravy = compute_gravy("RRRRR")
        assert gravy < -3.0

    def test_empty(self):
        assert compute_gravy("") == 0.0


class TestIsoelectricPoint:
    def test_neutral_peptide(self):
        pi = compute_isoelectric_point("AAAA")
        assert 5.0 < pi <= 7.5

    def test_basic_peptide(self):
        pi = compute_isoelectric_point("RRRR")
        assert pi > 10.0

    def test_acidic_peptide(self):
        pi = compute_isoelectric_point("DDDD")
        assert pi < 4.0


class TestInstabilityIndex:
    def test_returns_float(self):
        ii = compute_instability_index("AAAAAA")
        assert isinstance(ii, float)

    def test_short_sequence(self):
        assert compute_instability_index("A") == 0.0


class TestAromaticity:
    def test_no_aromatic(self):
        assert compute_aromaticity("AAAA") == 0.0

    def test_all_aromatic(self):
        assert compute_aromaticity("FWY") == pytest.approx(1.0, abs=0.01)

    def test_mixed(self):
        aro = compute_aromaticity("AAFAA")
        assert aro == pytest.approx(0.2, abs=0.01)


class TestProteaseSites:
    def test_trypsin(self):
        sites = count_protease_sites("AKRA")
        assert sites["trypsin"] == 2

    def test_proline_protection(self):
        # K before P should NOT be cleaved
        sites = count_protease_sites("AKPA")
        assert sites["trypsin"] == 0

    def test_chymotrypsin(self):
        sites = count_protease_sites("AFWA")
        assert sites["chymotrypsin"] == 2


class TestDrugLikeness:
    def test_good_peptide(self):
        score, violations = check_drug_likeness("ACDEFGHIKL")
        assert score >= 60
        assert len(violations) <= 2

    def test_too_short(self):
        score, violations = check_drug_likeness("ACE")
        assert any("length" in v for v in violations)

    def test_too_long(self):
        score, violations = check_drug_likeness("A" * 30)
        assert any("length" in v for v in violations)

    def test_hydrophobic_stretch(self):
        score, violations = check_drug_likeness("AVILLLLA")
        assert any("hydrophobic" in v for v in violations)


class TestPeptideProperties:
    def test_returns_all_keys(self):
        props = compute_peptide_properties("ACDEFGHIKL")
        expected_keys = [
            "sequence", "length", "molecular_weight", "net_charge",
            "gravy", "isoelectric_point", "instability_index",
            "aromaticity", "drug_likeness", "drug_likeness_violations",
            "protease_sites",
        ]
        for key in expected_keys:
            assert key in props


class TestSafeFloat:
    def test_valid(self):
        assert safe_float("3.14") == 3.14

    def test_invalid(self):
        assert safe_float("abc") == 0.0

    def test_none(self):
        assert safe_float(None, default=-1.0) == -1.0
