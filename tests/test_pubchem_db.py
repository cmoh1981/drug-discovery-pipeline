"""Tests for drugdiscovery.databases.pubchem_db module."""

from __future__ import annotations

from unittest.mock import MagicMock, call, patch

import pytest

from drugdiscovery.databases.pubchem_db import (
    _get_active_cids,
    _get_assay_ids,
    _get_compound_properties,
    _property_to_dict,
    search_by_target,
)


# ---------------------------------------------------------------------------
# Fixtures: reusable mock API responses
# ---------------------------------------------------------------------------

def _make_assay_response(aids: list[int]) -> dict:
    """Build a mock PubChem assay-IDs JSON response."""
    return {
        "InformationList": {
            "Information": [{"AID": aids}]
        }
    }


def _make_cids_response(cids: list[int]) -> dict:
    """Build a mock PubChem active-CIDs JSON response."""
    return {
        "InformationList": {
            "Information": [{"CID": cids}]
        }
    }


def _make_properties_response(properties: list[dict]) -> dict:
    """Build a mock PubChem compound properties JSON response."""
    return {
        "PropertyTable": {
            "Properties": properties,
        }
    }


def _mock_response(json_data: dict) -> MagicMock:
    """Create a mock requests.Response with .json() returning *json_data*."""
    resp = MagicMock()
    resp.json.return_value = json_data
    return resp


# ---------------------------------------------------------------------------
# Tests: search_by_target (main entry point)
# ---------------------------------------------------------------------------


class TestSearchByTarget:
    """Integration-level tests for the main search_by_target function."""

    @patch("drugdiscovery.databases.pubchem_db.time.sleep")
    @patch("drugdiscovery.databases.pubchem_db.get_with_retries")
    def test_full_workflow(self, mock_get: MagicMock, mock_sleep: MagicMock) -> None:
        """Successful search returns properly formatted compound dicts."""
        # 1st call: assay IDs
        aids_resp = _mock_response(_make_assay_response([1001, 1002]))
        # 2nd call: active CIDs for AID 1001
        cids_resp_1 = _mock_response(_make_cids_response([111, 222]))
        # 3rd call: active CIDs for AID 1002
        cids_resp_2 = _mock_response(_make_cids_response([333]))
        # 4th call: compound properties batch
        props_resp = _mock_response(_make_properties_response([
            {"CID": 111, "CanonicalSMILES": "CCO", "MolecularWeight": 46.07, "IUPACName": "ethanol"},
            {"CID": 222, "CanonicalSMILES": "CC(=O)O", "MolecularWeight": 60.05, "IUPACName": "acetic acid"},
            {"CID": 333, "CanonicalSMILES": "C(=O)O", "MolecularWeight": 46.03, "IUPACName": "formic acid"},
        ]))

        mock_get.side_effect = [aids_resp, cids_resp_1, cids_resp_2, props_resp]

        results = search_by_target("EGFR", limit=200)

        assert len(results) == 3
        assert results[0]["id"] == "CID111"
        assert results[0]["smiles"] == "CCO"
        assert results[0]["name"] == "ethanol"
        assert results[0]["source"] == "PubChem"
        assert results[0]["molecular_weight"] == 46.07
        assert results[0]["activity_type"] == "active"

        assert results[1]["id"] == "CID222"
        assert results[2]["id"] == "CID333"

    @patch("drugdiscovery.databases.pubchem_db.time.sleep")
    @patch("drugdiscovery.databases.pubchem_db.get_with_retries")
    def test_deduplicates_cids(self, mock_get: MagicMock, mock_sleep: MagicMock) -> None:
        """Duplicate CIDs across assays are deduplicated."""
        aids_resp = _mock_response(_make_assay_response([1001, 1002]))
        cids_resp_1 = _mock_response(_make_cids_response([111, 222]))
        cids_resp_2 = _mock_response(_make_cids_response([222, 333]))  # 222 is a duplicate
        props_resp = _mock_response(_make_properties_response([
            {"CID": 111, "CanonicalSMILES": "CCO", "MolecularWeight": 46.07, "IUPACName": "ethanol"},
            {"CID": 222, "CanonicalSMILES": "CC(=O)O", "MolecularWeight": 60.05, "IUPACName": "acetic acid"},
            {"CID": 333, "CanonicalSMILES": "C(=O)O", "MolecularWeight": 46.03, "IUPACName": "formic acid"},
        ]))

        mock_get.side_effect = [aids_resp, cids_resp_1, cids_resp_2, props_resp]

        results = search_by_target("BRAF", limit=200)

        assert len(results) == 3
        # Verify the property fetch used deduplicated CIDs (111, 222, 333)
        property_call = mock_get.call_args_list[3]
        url_arg = property_call[0][0]
        assert "111,222,333" in url_arg

    @patch("drugdiscovery.databases.pubchem_db.time.sleep")
    @patch("drugdiscovery.databases.pubchem_db.get_with_retries")
    def test_respects_limit(self, mock_get: MagicMock, mock_sleep: MagicMock) -> None:
        """Result count is capped at the limit parameter."""
        aids_resp = _mock_response(_make_assay_response([1001]))
        cids_resp = _mock_response(_make_cids_response([111, 222, 333, 444, 555]))
        props_resp = _mock_response(_make_properties_response([
            {"CID": 111, "CanonicalSMILES": "CCO", "MolecularWeight": 46.07, "IUPACName": "ethanol"},
            {"CID": 222, "CanonicalSMILES": "CC(=O)O", "MolecularWeight": 60.05, "IUPACName": "acetic acid"},
        ]))

        mock_get.side_effect = [aids_resp, cids_resp, props_resp]

        results = search_by_target("TP53", limit=2)

        # The property fetch URL should only contain 2 CIDs
        property_call = mock_get.call_args_list[2]
        url_arg = property_call[0][0]
        assert "111,222" in url_arg
        assert "333" not in url_arg


# ---------------------------------------------------------------------------
# Tests: empty results handling
# ---------------------------------------------------------------------------


class TestEmptyResults:
    """Verify graceful empty-list returns for various edge cases."""

    @patch("drugdiscovery.databases.pubchem_db.get_with_retries")
    def test_no_assays_found(self, mock_get: MagicMock) -> None:
        """Return empty list when gene has no associated assays."""
        mock_get.return_value = _mock_response(_make_assay_response([]))

        results = search_by_target("NONEXISTENT_GENE")
        assert results == []

    @patch("drugdiscovery.databases.pubchem_db.time.sleep")
    @patch("drugdiscovery.databases.pubchem_db.get_with_retries")
    def test_no_active_compounds(self, mock_get: MagicMock, mock_sleep: MagicMock) -> None:
        """Return empty list when assays exist but no active CIDs found."""
        aids_resp = _mock_response(_make_assay_response([1001, 1002]))
        cids_resp_empty = _mock_response(_make_cids_response([]))

        mock_get.side_effect = [aids_resp, cids_resp_empty, cids_resp_empty]

        results = search_by_target("RARE_TARGET")
        assert results == []

    @patch("drugdiscovery.databases.pubchem_db.get_with_retries")
    def test_empty_information_list(self, mock_get: MagicMock) -> None:
        """Handle edge case where InformationList has empty Information."""
        mock_get.return_value = _mock_response({
            "InformationList": {"Information": [{}]}
        })

        results = search_by_target("EDGE_CASE")
        assert results == []


# ---------------------------------------------------------------------------
# Tests: rate limiting
# ---------------------------------------------------------------------------


class TestRateLimiting:
    """Verify rate-limiting delays between API calls."""

    @patch("drugdiscovery.databases.pubchem_db.time.sleep")
    @patch("drugdiscovery.databases.pubchem_db.get_with_retries")
    def test_sleep_between_assay_cid_fetches(
        self, mock_get: MagicMock, mock_sleep: MagicMock
    ) -> None:
        """A sleep(0.25) is called before each active-CID fetch."""
        aids_resp = _mock_response(_make_assay_response([1001, 1002, 1003]))
        cids_resp = _mock_response(_make_cids_response([111]))
        props_resp = _mock_response(_make_properties_response([
            {"CID": 111, "CanonicalSMILES": "CCO", "MolecularWeight": 46.07, "IUPACName": "ethanol"},
        ]))

        mock_get.side_effect = [aids_resp, cids_resp, cids_resp, cids_resp, props_resp]

        search_by_target("KRAS")

        # 3 sleeps for CID fetches + 1 sleep for property batch = 4 total
        sleep_calls = [c for c in mock_sleep.call_args_list if c == call(0.25)]
        assert len(sleep_calls) >= 3

    @patch("drugdiscovery.databases.pubchem_db.time.sleep")
    @patch("drugdiscovery.databases.pubchem_db.get_with_retries")
    def test_sleep_between_property_batches(
        self, mock_get: MagicMock, mock_sleep: MagicMock
    ) -> None:
        """A sleep(0.25) is called between compound property batch fetches."""
        # Use 4 assays each returning 50 CIDs = 200 total, to trigger
        # 2 property batches (batch_size=100). _get_active_cids caps at
        # max_cids=50, so we need multiple assays to exceed 100 CIDs.
        cids_a = list(range(1, 51))
        cids_b = list(range(51, 101))
        cids_c = list(range(101, 151))
        cids_d = list(range(151, 201))

        aids_resp = _mock_response(_make_assay_response([1001, 1002, 1003, 1004]))
        cids_resp_a = _mock_response(_make_cids_response(cids_a))
        cids_resp_b = _mock_response(_make_cids_response(cids_b))
        cids_resp_c = _mock_response(_make_cids_response(cids_c))
        cids_resp_d = _mock_response(_make_cids_response(cids_d))

        all_cids = cids_a + cids_b + cids_c + cids_d
        props_batch_1 = _make_properties_response([
            {"CID": c, "CanonicalSMILES": "C", "MolecularWeight": 12.0, "IUPACName": "methane"}
            for c in all_cids[:100]
        ])
        props_batch_2 = _make_properties_response([
            {"CID": c, "CanonicalSMILES": "C", "MolecularWeight": 12.0, "IUPACName": "methane"}
            for c in all_cids[100:]
        ])

        mock_get.side_effect = [
            aids_resp,
            cids_resp_a, cids_resp_b, cids_resp_c, cids_resp_d,
            _mock_response(props_batch_1),
            _mock_response(props_batch_2),
        ]

        results = search_by_target("LARGE_TARGET", limit=200)

        assert len(results) == 200
        # At least 2 sleeps for property batches + 4 for CID fetches
        sleep_calls = [c for c in mock_sleep.call_args_list if c == call(0.25)]
        assert len(sleep_calls) >= 6


# ---------------------------------------------------------------------------
# Tests: error handling (API failures)
# ---------------------------------------------------------------------------


class TestErrorHandling:
    """Verify graceful degradation when API calls fail."""

    @patch("drugdiscovery.databases.pubchem_db.get_with_retries")
    def test_assay_api_failure_returns_empty(self, mock_get: MagicMock) -> None:
        """Return empty list when the assay lookup request raises."""
        mock_get.side_effect = Exception("Connection refused")

        results = search_by_target("EGFR")
        assert results == []

    @patch("drugdiscovery.databases.pubchem_db.time.sleep")
    @patch("drugdiscovery.databases.pubchem_db.get_with_retries")
    def test_cid_api_failure_continues(
        self, mock_get: MagicMock, mock_sleep: MagicMock
    ) -> None:
        """When one assay's CID fetch fails, other assays still proceed."""
        aids_resp = _mock_response(_make_assay_response([1001, 1002]))
        # AID 1001 fails
        cids_error = Exception("Timeout")
        # AID 1002 succeeds
        cids_resp_ok = _mock_response(_make_cids_response([444]))
        props_resp = _mock_response(_make_properties_response([
            {"CID": 444, "CanonicalSMILES": "CC", "MolecularWeight": 30.07, "IUPACName": "ethane"},
        ]))

        mock_get.side_effect = [aids_resp, cids_error, cids_resp_ok, props_resp]

        results = search_by_target("MIXED_TARGET")

        assert len(results) == 1
        assert results[0]["id"] == "CID444"

    @patch("drugdiscovery.databases.pubchem_db.time.sleep")
    @patch("drugdiscovery.databases.pubchem_db.get_with_retries")
    def test_property_batch_failure_continues(
        self, mock_get: MagicMock, mock_sleep: MagicMock
    ) -> None:
        """When a property batch fails, other batches still return results."""
        # Use 4 assays each returning 50 CIDs = 200 total, triggering
        # 2 property batches. _get_active_cids caps at max_cids=50.
        cids_a = list(range(1, 51))
        cids_b = list(range(51, 101))
        cids_c = list(range(101, 151))
        cids_d = list(range(151, 201))

        aids_resp = _mock_response(_make_assay_response([1001, 1002, 1003, 1004]))
        cids_resp_a = _mock_response(_make_cids_response(cids_a))
        cids_resp_b = _mock_response(_make_cids_response(cids_b))
        cids_resp_c = _mock_response(_make_cids_response(cids_c))
        cids_resp_d = _mock_response(_make_cids_response(cids_d))

        # First property batch (CIDs 1-100) fails, second (CIDs 101-200) succeeds
        props_batch_2 = _make_properties_response([
            {"CID": c, "CanonicalSMILES": "C", "MolecularWeight": 12.0, "IUPACName": "methane"}
            for c in range(101, 201)
        ])

        mock_get.side_effect = [
            aids_resp,
            cids_resp_a, cids_resp_b, cids_resp_c, cids_resp_d,
            Exception("Server error"),
            _mock_response(props_batch_2),
        ]

        results = search_by_target("PARTIAL_TARGET", limit=200)

        # Only second batch (100 compounds) should succeed
        assert len(results) == 100

    @patch("drugdiscovery.databases.pubchem_db.get_with_retries")
    def test_malformed_json_returns_empty(self, mock_get: MagicMock) -> None:
        """Return empty list on malformed JSON response."""
        mock_get.return_value = _mock_response({"unexpected": "format"})

        results = search_by_target("MALFORMED")
        assert results == []


# ---------------------------------------------------------------------------
# Tests: internal helpers
# ---------------------------------------------------------------------------


class TestGetAssayIds:
    """Unit tests for _get_assay_ids."""

    @patch("drugdiscovery.databases.pubchem_db.get_with_retries")
    def test_returns_limited_aids(self, mock_get: MagicMock) -> None:
        """Respect max_assays parameter."""
        aids = list(range(1, 25))
        mock_get.return_value = _mock_response(_make_assay_response(aids))

        result = _get_assay_ids("EGFR", max_assays=5)
        assert len(result) == 5
        assert result == [1, 2, 3, 4, 5]

    @patch("drugdiscovery.databases.pubchem_db.get_with_retries")
    def test_failure_returns_empty(self, mock_get: MagicMock) -> None:
        """Return empty list on exception."""
        mock_get.side_effect = Exception("Network error")
        assert _get_assay_ids("BRAF") == []


class TestGetActiveCids:
    """Unit tests for _get_active_cids."""

    @patch("drugdiscovery.databases.pubchem_db.get_with_retries")
    def test_returns_limited_cids(self, mock_get: MagicMock) -> None:
        """Respect max_cids parameter."""
        cids = list(range(1, 200))
        mock_get.return_value = _mock_response(_make_cids_response(cids))

        result = _get_active_cids(12345, max_cids=10)
        assert len(result) == 10

    @patch("drugdiscovery.databases.pubchem_db.get_with_retries")
    def test_passes_cids_type_param(self, mock_get: MagicMock) -> None:
        """Verify cids_type=active is passed as a query parameter."""
        mock_get.return_value = _mock_response(_make_cids_response([111]))

        _get_active_cids(99999)

        _, kwargs = mock_get.call_args
        assert kwargs["params"] == {"cids_type": "active"}

    @patch("drugdiscovery.databases.pubchem_db.get_with_retries")
    def test_failure_returns_empty(self, mock_get: MagicMock) -> None:
        """Return empty list on exception."""
        mock_get.side_effect = Exception("502 Bad Gateway")
        assert _get_active_cids(12345) == []


class TestGetCompoundProperties:
    """Unit tests for _get_compound_properties."""

    @patch("drugdiscovery.databases.pubchem_db.time.sleep")
    @patch("drugdiscovery.databases.pubchem_db.get_with_retries")
    def test_single_batch(self, mock_get: MagicMock, mock_sleep: MagicMock) -> None:
        """Fewer than 100 CIDs are fetched in a single batch."""
        mock_get.return_value = _mock_response(_make_properties_response([
            {"CID": 111, "CanonicalSMILES": "CCO", "MolecularWeight": 46.07, "IUPACName": "ethanol"},
        ]))

        results = _get_compound_properties([111])
        assert len(results) == 1
        assert mock_get.call_count == 1

    @patch("drugdiscovery.databases.pubchem_db.time.sleep")
    @patch("drugdiscovery.databases.pubchem_db.get_with_retries")
    def test_multiple_batches(self, mock_get: MagicMock, mock_sleep: MagicMock) -> None:
        """More than 100 CIDs are split into multiple batches."""
        cids = list(range(1, 251))  # 250 CIDs -> 3 batches

        batch_resp = _mock_response(_make_properties_response([
            {"CID": c, "CanonicalSMILES": "C", "MolecularWeight": 12.0, "IUPACName": "methane"}
            for c in range(1, 101)
        ]))
        mock_get.return_value = batch_resp

        results = _get_compound_properties(cids)

        assert mock_get.call_count == 3


class TestPropertyToDict:
    """Unit tests for _property_to_dict."""

    def test_complete_record(self) -> None:
        """All fields are correctly mapped."""
        prop = {
            "CID": 12345,
            "CanonicalSMILES": "c1ccccc1",
            "MolecularWeight": 78.11,
            "IUPACName": "benzene",
        }
        result = _property_to_dict(prop)

        assert result == {
            "id": "CID12345",
            "smiles": "c1ccccc1",
            "name": "benzene",
            "source": "PubChem",
            "molecular_weight": 78.11,
            "activity_type": "active",
            "activity_value": None,
        }

    def test_missing_fields_use_defaults(self) -> None:
        """Missing optional fields fall back to safe defaults."""
        result = _property_to_dict({})

        assert result["id"] == ""
        assert result["smiles"] == ""
        assert result["name"] == ""
        assert result["source"] == "PubChem"
        assert result["molecular_weight"] is None
