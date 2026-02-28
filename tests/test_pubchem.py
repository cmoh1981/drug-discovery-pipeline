"""Tests for drugdiscovery.tools.pubchem module."""

from unittest.mock import MagicMock, patch

from drugdiscovery.tools.pubchem import get_iupac_name, resolve_iupac_names
from drugdiscovery.types import Candidate


class TestGetIupacName:
    @patch("drugdiscovery.tools.pubchem.get_with_retries")
    def test_valid_response(self, mock_get):
        mock_resp = MagicMock()
        mock_resp.json.return_value = {
            "PropertyTable": {
                "Properties": [{"CID": 12345, "IUPACName": "pyridine-3-sulfonamide"}]
            }
        }
        mock_get.return_value = mock_resp

        result = get_iupac_name("NS(=O)(=O)c1cccnc1")
        assert result == "pyridine-3-sulfonamide"
        mock_get.assert_called_once()

    @patch("drugdiscovery.tools.pubchem.get_with_retries")
    def test_api_error_returns_empty(self, mock_get):
        mock_get.side_effect = Exception("API down")
        result = get_iupac_name("CC(=O)O")
        assert result == ""

    def test_empty_smiles_returns_empty(self):
        assert get_iupac_name("") == ""


class TestResolveIupacNames:
    @patch("drugdiscovery.tools.pubchem.get_iupac_name")
    @patch("drugdiscovery.tools.pubchem.time.sleep")
    def test_skips_peptide_candidates(self, mock_sleep, mock_lookup):
        candidates = [
            Candidate(candidate_id="PEP_001", modality="peptide", sequence="ACDEFG"),
            Candidate(candidate_id="PEP_002", modality="peptide", sequence="GHIKLM"),
        ]
        resolve_iupac_names(candidates)
        mock_lookup.assert_not_called()
        assert candidates[0].iupac_name == ""

    @patch("drugdiscovery.tools.pubchem.get_iupac_name")
    @patch("drugdiscovery.tools.pubchem.time.sleep")
    def test_caches_duplicate_smiles(self, mock_sleep, mock_lookup):
        mock_lookup.return_value = "acetic acid"
        candidates = [
            Candidate(candidate_id="SM_001", modality="small_molecule", smiles="CC(=O)O"),
            Candidate(candidate_id="SM_002", modality="small_molecule", smiles="CC(=O)O",
                      parent_id="SM_001", modification="methylation"),
        ]
        resolve_iupac_names(candidates)
        # Only one actual API call despite two candidates with same SMILES
        mock_lookup.assert_called_once_with("CC(=O)O")
        assert candidates[0].iupac_name == "acetic acid"
        assert candidates[1].iupac_name == "acetic acid"

    @patch("drugdiscovery.tools.pubchem.get_iupac_name")
    @patch("drugdiscovery.tools.pubchem.time.sleep")
    def test_resolves_small_molecules(self, mock_sleep, mock_lookup):
        mock_lookup.return_value = "ethanol"
        candidates = [
            Candidate(candidate_id="SM_001", modality="small_molecule", smiles="CCO"),
        ]
        resolve_iupac_names(candidates)
        assert candidates[0].iupac_name == "ethanol"


class TestCandidateToDictIncludesIupac:
    def test_iupac_in_dict(self):
        c = Candidate(
            candidate_id="SM_001",
            modality="small_molecule",
            smiles="CCO",
            iupac_name="ethanol",
        )
        d = c.to_dict()
        assert "iupac_name" in d
        assert d["iupac_name"] == "ethanol"

    def test_iupac_default_empty_in_dict(self):
        c = Candidate(candidate_id="SM_002")
        d = c.to_dict()
        assert d["iupac_name"] == ""
