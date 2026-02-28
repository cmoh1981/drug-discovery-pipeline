"""Tests for drugdiscovery.databases.bindingdb module."""

from unittest.mock import MagicMock, patch

from drugdiscovery.databases.bindingdb import (
    _best_affinity,
    _parse_binding_response,
    _safe_affinity,
    search_by_target,
    search_by_uniprot,
)


# ---------------------------------------------------------------------------
# _safe_affinity
# ---------------------------------------------------------------------------


class TestSafeAffinity:
    def test_numeric_float(self):
        assert _safe_affinity(10.5) == 10.5

    def test_numeric_int(self):
        assert _safe_affinity(100) == 100.0

    def test_string_number(self):
        assert _safe_affinity("42.3") == 42.3

    def test_greater_than_prefix(self):
        assert _safe_affinity(">10000") == 10000.0

    def test_less_than_prefix(self):
        assert _safe_affinity("<1") == 1.0

    def test_gte_prefix(self):
        assert _safe_affinity(">=500") == 500.0

    def test_lte_prefix(self):
        assert _safe_affinity("<=250") == 250.0

    def test_tilde_prefix(self):
        assert _safe_affinity("~100") == 100.0

    def test_none_returns_none(self):
        assert _safe_affinity(None) is None

    def test_empty_string_returns_none(self):
        assert _safe_affinity("") is None

    def test_whitespace_returns_none(self):
        assert _safe_affinity("   ") is None

    def test_non_numeric_returns_none(self):
        assert _safe_affinity("not-a-number") is None

    def test_negative_returns_none(self):
        assert _safe_affinity(-5.0) is None


# ---------------------------------------------------------------------------
# _best_affinity
# ---------------------------------------------------------------------------


class TestBestAffinity:
    def test_ki_only(self):
        assert _best_affinity(10.0, None, None, None) == ("Ki", 10.0)

    def test_kd_only(self):
        assert _best_affinity(None, 5.0, None, None) == ("Kd", 5.0)

    def test_ic50_only(self):
        assert _best_affinity(None, None, 100.0, None) == ("IC50", 100.0)

    def test_ec50_only(self):
        assert _best_affinity(None, None, None, 50.0) == ("EC50", 50.0)

    def test_lowest_wins(self):
        assert _best_affinity(100.0, 5.0, 50.0, 200.0) == ("Kd", 5.0)

    def test_all_none_returns_empty(self):
        assert _best_affinity(None, None, None, None) == ("", None)

    def test_tie_goes_to_first(self):
        # Ki and Kd both 10.0 — Ki comes first in candidate list
        atype, aval = _best_affinity(10.0, 10.0, None, None)
        assert aval == 10.0
        assert atype in ("Ki", "Kd")


# ---------------------------------------------------------------------------
# _parse_binding_response
# ---------------------------------------------------------------------------


class TestParseBindingResponse:
    def test_uniprot_envelope(self):
        data = {
            "getLigandsByUniprotsResponse": {
                "affinities": [
                    {
                        "monomerid": "BDB_12345",
                        "smiles": "CCO",
                        "ligand_name": "ethanol",
                        "ki": "100",
                        "kd": None,
                        "ic50": "200",
                        "ec50": None,
                    }
                ]
            }
        }
        results = _parse_binding_response(data, limit=10)
        assert len(results) == 1
        assert results[0]["id"] == "BDB_12345"
        assert results[0]["smiles"] == "CCO"
        assert results[0]["source"] == "BindingDB"
        assert results[0]["ki"] == 100.0
        assert results[0]["ic50"] == 200.0
        assert results[0]["affinity_type"] == "Ki"
        assert results[0]["affinity_value"] == 100.0

    def test_target_envelope(self):
        data = {
            "getLigandsByTargetResponse": {
                "affinities": [
                    {
                        "monomerid": "BDB_999",
                        "smiles": "c1ccccc1",
                        "ligand_name": "benzene",
                        "ki": None,
                        "kd": "50",
                        "ic50": None,
                        "ec50": None,
                    }
                ]
            }
        }
        results = _parse_binding_response(data, limit=10)
        assert len(results) == 1
        assert results[0]["affinity_type"] == "Kd"
        assert results[0]["affinity_value"] == 50.0

    def test_bare_list(self):
        data = [
            {
                "monomerid": "1",
                "smiles": "C",
                "ligand_name": "methane",
                "ki": "10",
                "kd": None,
                "ic50": None,
                "ec50": None,
            },
            {
                "monomerid": "2",
                "smiles": "CC",
                "ligand_name": "ethane",
                "ki": None,
                "kd": None,
                "ic50": "500",
                "ec50": None,
            },
        ]
        results = _parse_binding_response(data, limit=10)
        assert len(results) == 2

    def test_single_record_as_dict(self):
        """BindingDB may return a single record as a dict instead of list."""
        data = {
            "getLigandsByUniprotsResponse": {
                "affinities": {
                    "monomerid": "SINGLE",
                    "smiles": "O",
                    "ligand_name": "water",
                    "ki": "1",
                    "kd": None,
                    "ic50": None,
                    "ec50": None,
                }
            }
        }
        results = _parse_binding_response(data, limit=10)
        assert len(results) == 1
        assert results[0]["id"] == "SINGLE"

    def test_limit_truncates(self):
        records = [
            {"monomerid": str(i), "smiles": "C", "ligand_name": f"mol_{i}"}
            for i in range(50)
        ]
        results = _parse_binding_response(records, limit=5)
        assert len(results) == 5

    def test_empty_response(self):
        assert _parse_binding_response({}, limit=10) == []
        assert _parse_binding_response([], limit=10) == []

    def test_no_affinity_data(self):
        data = [{"monomerid": "X", "smiles": "C", "ligand_name": "test"}]
        results = _parse_binding_response(data, limit=10)
        assert len(results) == 1
        assert results[0]["affinity_type"] == ""
        assert results[0]["affinity_value"] is None


# ---------------------------------------------------------------------------
# search_by_uniprot
# ---------------------------------------------------------------------------


class TestSearchByUniprot:
    @patch("drugdiscovery.databases.bindingdb.get_with_retries")
    def test_successful_search(self, mock_get):
        mock_resp = MagicMock()
        mock_resp.json.return_value = {
            "getLigandsByUniprotsResponse": {
                "affinities": [
                    {
                        "monomerid": "BDB_001",
                        "smiles": "CCO",
                        "ligand_name": "ethanol",
                        "ki": "50",
                        "kd": None,
                        "ic50": None,
                        "ec50": None,
                    }
                ]
            }
        }
        mock_get.return_value = mock_resp

        results = search_by_uniprot("P00533")
        assert len(results) == 1
        assert results[0]["id"] == "BDB_001"
        assert results[0]["source"] == "BindingDB"
        mock_get.assert_called_once()
        call_kwargs = mock_get.call_args
        assert call_kwargs[1]["params"]["uniprot"] == "P00533"
        assert call_kwargs[1]["params"]["response"] == "json"
        assert call_kwargs[1]["timeout"] == 60

    @patch("drugdiscovery.databases.bindingdb.get_with_retries")
    def test_api_error_returns_empty(self, mock_get):
        mock_get.side_effect = Exception("Connection refused")
        results = search_by_uniprot("P00533")
        assert results == []

    def test_empty_uniprot_returns_empty(self):
        assert search_by_uniprot("") == []

    @patch("drugdiscovery.databases.bindingdb.get_with_retries")
    def test_limit_is_applied(self, mock_get):
        records = [
            {
                "monomerid": f"BDB_{i}",
                "smiles": "C",
                "ligand_name": f"mol_{i}",
                "ki": str(i * 10),
            }
            for i in range(20)
        ]
        mock_resp = MagicMock()
        mock_resp.json.return_value = {
            "getLigandsByUniprotsResponse": {"affinities": records}
        }
        mock_get.return_value = mock_resp

        results = search_by_uniprot("P00533", limit=5)
        assert len(results) == 5


# ---------------------------------------------------------------------------
# search_by_target
# ---------------------------------------------------------------------------


class TestSearchByTarget:
    @patch("drugdiscovery.databases.bindingdb.get_with_retries")
    def test_successful_search(self, mock_get):
        mock_resp = MagicMock()
        mock_resp.json.return_value = {
            "getLigandsByTargetResponse": {
                "affinities": [
                    {
                        "monomerid": "BDB_100",
                        "smiles": "c1ccccc1",
                        "ligand_name": "benzene",
                        "ki": None,
                        "kd": "25",
                        "ic50": "100",
                        "ec50": None,
                    }
                ]
            }
        }
        mock_get.return_value = mock_resp

        results = search_by_target("EGFR")
        assert len(results) == 1
        assert results[0]["id"] == "BDB_100"
        assert results[0]["kd"] == 25.0
        assert results[0]["ic50"] == 100.0
        assert results[0]["affinity_type"] == "Kd"
        assert results[0]["affinity_value"] == 25.0
        mock_get.assert_called_once()
        call_kwargs = mock_get.call_args
        assert call_kwargs[1]["params"]["target"] == "EGFR"
        assert call_kwargs[1]["params"]["response"] == "json"

    @patch("drugdiscovery.databases.bindingdb.get_with_retries")
    def test_api_error_returns_empty(self, mock_get):
        mock_get.side_effect = Exception("Timeout")
        results = search_by_target("EGFR")
        assert results == []

    def test_empty_target_returns_empty(self):
        assert search_by_target("") == []

    @patch("drugdiscovery.databases.bindingdb.get_with_retries")
    def test_multiple_results_with_best_affinity(self, mock_get):
        mock_resp = MagicMock()
        mock_resp.json.return_value = {
            "getLigandsByTargetResponse": {
                "affinities": [
                    {
                        "monomerid": "BDB_A",
                        "smiles": "CCO",
                        "ligand_name": "compound_a",
                        "ki": ">10000",
                        "kd": None,
                        "ic50": "500",
                        "ec50": None,
                    },
                    {
                        "monomerid": "BDB_B",
                        "smiles": "CCN",
                        "ligand_name": "compound_b",
                        "ki": "1.5",
                        "kd": "3.0",
                        "ic50": None,
                        "ec50": None,
                    },
                ]
            }
        }
        mock_get.return_value = mock_resp

        results = search_by_target("BRAF")
        assert len(results) == 2

        # First compound: IC50 (500) is best since Ki is >10000
        assert results[0]["ki"] == 10000.0
        assert results[0]["ic50"] == 500.0
        assert results[0]["affinity_type"] == "IC50"
        assert results[0]["affinity_value"] == 500.0

        # Second compound: Ki (1.5) is strongest
        assert results[1]["ki"] == 1.5
        assert results[1]["kd"] == 3.0
        assert results[1]["affinity_type"] == "Ki"
        assert results[1]["affinity_value"] == 1.5
