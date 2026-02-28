"""Tests for drugdiscovery.databases.opentargets module.

Tests cover:
  - Target ID resolution from gene symbol
  - Target-disease association queries
  - Known drugs retrieval
  - Disease gene aggregation
  - Error handling and edge cases
"""

from unittest.mock import MagicMock, patch, call

from drugdiscovery.databases.opentargets import (
    get_target_diseases,
    get_known_drugs,
    get_disease_genes,
    _resolve_target_id,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _mock_graphql_response(data: dict) -> MagicMock:
    """Create a mock response object returning the given data as JSON."""
    mock_resp = MagicMock()
    mock_resp.json.return_value = data
    return mock_resp


# ---------------------------------------------------------------------------
# Target ID Resolution Tests
# ---------------------------------------------------------------------------

class TestResolveTargetId:
    """Test _resolve_target_id mapping gene symbols to Ensembl IDs."""

    @patch("drugdiscovery.databases.opentargets.post_with_retries")
    def test_resolves_first_ensg_hit(self, mock_post):
        mock_post.return_value = _mock_graphql_response({
            "data": {
                "search": {
                    "hits": [
                        {"id": "ENSG00000146648"},
                        {"id": "ENSG00000099999"},
                    ]
                }
            }
        })

        result = _resolve_target_id("EGFR")
        assert result == "ENSG00000146648"
        mock_post.assert_called_once()

    @patch("drugdiscovery.databases.opentargets.post_with_retries")
    def test_skips_non_ensg_first_hit(self, mock_post):
        mock_post.return_value = _mock_graphql_response({
            "data": {
                "search": {
                    "hits": [
                        {"id": "EFO_0000311"},  # disease ID, not target
                        {"id": "ENSG00000146648"},
                    ]
                }
            }
        })

        result = _resolve_target_id("EGFR")
        assert result == "ENSG00000146648"

    @patch("drugdiscovery.databases.opentargets.post_with_retries")
    def test_returns_none_when_no_hits(self, mock_post):
        mock_post.return_value = _mock_graphql_response({
            "data": {
                "search": {
                    "hits": []
                }
            }
        })

        result = _resolve_target_id("NONEXISTENT_GENE")
        assert result is None

    @patch("drugdiscovery.databases.opentargets.post_with_retries")
    def test_returns_none_when_no_ensg_ids(self, mock_post):
        mock_post.return_value = _mock_graphql_response({
            "data": {
                "search": {
                    "hits": [
                        {"id": "EFO_0000311"},
                        {"id": "CHEMBL203"},
                    ]
                }
            }
        })

        result = _resolve_target_id("AMBIGUOUS")
        assert result is None


# ---------------------------------------------------------------------------
# Target-Disease Association Tests
# ---------------------------------------------------------------------------

class TestGetTargetDiseases:
    """Test get_target_diseases fetching disease associations."""

    @patch("drugdiscovery.databases.opentargets.post_with_retries")
    def test_returns_disease_associations(self, mock_post):
        # First call: target search; second call: associated diseases
        mock_post.side_effect = [
            _mock_graphql_response({
                "data": {
                    "search": {
                        "hits": [{"id": "ENSG00000146648"}]
                    }
                }
            }),
            _mock_graphql_response({
                "data": {
                    "target": {
                        "id": "ENSG00000146648",
                        "approvedSymbol": "EGFR",
                        "associatedDiseases": {
                            "rows": [
                                {
                                    "disease": {
                                        "id": "EFO_0000311",
                                        "name": "lung carcinoma",
                                        "therapeuticAreas": [
                                            {"id": "OTAR_0000018", "name": "neoplasm"}
                                        ],
                                    },
                                    "score": 0.85,
                                },
                                {
                                    "disease": {
                                        "id": "EFO_0000616",
                                        "name": "glioblastoma",
                                        "therapeuticAreas": [
                                            {"id": "OTAR_0000018", "name": "neoplasm"}
                                        ],
                                    },
                                    "score": 0.72,
                                },
                            ]
                        },
                    }
                }
            }),
        ]

        results = get_target_diseases("EGFR", limit=5)
        assert len(results) == 2
        assert results[0]["disease_id"] == "EFO_0000311"
        assert results[0]["disease_name"] == "lung carcinoma"
        assert results[0]["score"] == 0.85
        assert "neoplasm" in results[0]["therapeutic_areas"]
        assert results[1]["disease_name"] == "glioblastoma"

    @patch("drugdiscovery.databases.opentargets.post_with_retries")
    def test_returns_empty_on_unknown_gene(self, mock_post):
        mock_post.return_value = _mock_graphql_response({
            "data": {
                "search": {
                    "hits": []
                }
            }
        })

        results = get_target_diseases("NONEXISTENT_GENE_XYZ")
        assert results == []

    @patch("drugdiscovery.databases.opentargets.post_with_retries")
    def test_returns_empty_on_api_error(self, mock_post):
        mock_post.side_effect = Exception("API timeout")

        results = get_target_diseases("EGFR")
        assert results == []

    @patch("drugdiscovery.databases.opentargets.post_with_retries")
    def test_handles_missing_therapeutic_areas(self, mock_post):
        mock_post.side_effect = [
            _mock_graphql_response({
                "data": {
                    "search": {
                        "hits": [{"id": "ENSG00000146648"}]
                    }
                }
            }),
            _mock_graphql_response({
                "data": {
                    "target": {
                        "id": "ENSG00000146648",
                        "approvedSymbol": "EGFR",
                        "associatedDiseases": {
                            "rows": [
                                {
                                    "disease": {
                                        "id": "EFO_0000311",
                                        "name": "lung carcinoma",
                                    },
                                    "score": 0.85,
                                },
                            ]
                        },
                    }
                }
            }),
        ]

        results = get_target_diseases("EGFR")
        assert len(results) == 1
        assert results[0]["therapeutic_areas"] == []


# ---------------------------------------------------------------------------
# Known Drugs Tests
# ---------------------------------------------------------------------------

class TestGetKnownDrugs:
    """Test get_known_drugs fetching drug/clinical compound data."""

    @patch("drugdiscovery.databases.opentargets.post_with_retries")
    def test_returns_known_drugs(self, mock_post):
        mock_post.side_effect = [
            _mock_graphql_response({
                "data": {
                    "search": {
                        "hits": [{"id": "ENSG00000146648"}]
                    }
                }
            }),
            _mock_graphql_response({
                "data": {
                    "target": {
                        "knownDrugs": {
                            "rows": [
                                {
                                    "drugId": "CHEMBL553",
                                    "prefName": "erlotinib",
                                    "drugType": "Small molecule",
                                    "mechanismOfAction": "EGFR inhibitor",
                                    "phase": 4,
                                },
                                {
                                    "drugId": "CHEMBL1201583",
                                    "prefName": "gefitinib",
                                    "drugType": "Small molecule",
                                    "mechanismOfAction": "EGFR inhibitor",
                                    "phase": 4,
                                },
                            ]
                        }
                    }
                }
            }),
        ]

        results = get_known_drugs("EGFR", limit=10)
        assert len(results) == 2
        assert results[0]["id"] == "CHEMBL553"
        assert results[0]["drug_name"] == "erlotinib"
        assert results[0]["drug_type"] == "Small molecule"
        assert results[0]["phase"] == 4
        assert results[0]["mechanism_of_action"] == "EGFR inhibitor"
        assert results[0]["source"] == "OpenTargets"

    @patch("drugdiscovery.databases.opentargets.post_with_retries")
    def test_deduplicates_drug_ids(self, mock_post):
        mock_post.side_effect = [
            _mock_graphql_response({
                "data": {
                    "search": {
                        "hits": [{"id": "ENSG00000146648"}]
                    }
                }
            }),
            _mock_graphql_response({
                "data": {
                    "target": {
                        "knownDrugs": {
                            "rows": [
                                {
                                    "drugId": "CHEMBL553",
                                    "prefName": "erlotinib",
                                    "drugType": "Small molecule",
                                    "mechanismOfAction": "EGFR inhibitor",
                                    "phase": 4,
                                },
                                {
                                    "drugId": "CHEMBL553",
                                    "prefName": "erlotinib",
                                    "drugType": "Small molecule",
                                    "mechanismOfAction": "EGFR blocker",
                                    "phase": 3,
                                },
                            ]
                        }
                    }
                }
            }),
        ]

        results = get_known_drugs("EGFR")
        assert len(results) == 1
        assert results[0]["drug_name"] == "erlotinib"

    @patch("drugdiscovery.databases.opentargets.post_with_retries")
    def test_returns_empty_on_unknown_gene(self, mock_post):
        mock_post.return_value = _mock_graphql_response({
            "data": {
                "search": {
                    "hits": []
                }
            }
        })

        results = get_known_drugs("NONEXISTENT_GENE_XYZ")
        assert results == []

    @patch("drugdiscovery.databases.opentargets.post_with_retries")
    def test_returns_empty_on_api_error(self, mock_post):
        mock_post.side_effect = Exception("Connection refused")

        results = get_known_drugs("EGFR")
        assert results == []

    @patch("drugdiscovery.databases.opentargets.post_with_retries")
    def test_uses_drug_id_as_name_fallback(self, mock_post):
        mock_post.side_effect = [
            _mock_graphql_response({
                "data": {
                    "search": {
                        "hits": [{"id": "ENSG00000146648"}]
                    }
                }
            }),
            _mock_graphql_response({
                "data": {
                    "target": {
                        "knownDrugs": {
                            "rows": [
                                {
                                    "drugId": "CHEMBL999999",
                                    "prefName": None,
                                    "drugType": "Antibody",
                                    "mechanismOfAction": "",
                                    "phase": 1,
                                },
                            ]
                        }
                    }
                }
            }),
        ]

        results = get_known_drugs("EGFR")
        assert len(results) == 1
        assert results[0]["drug_name"] == "CHEMBL999999"


# ---------------------------------------------------------------------------
# Disease Genes Tests
# ---------------------------------------------------------------------------

class TestGetDiseaseGenes:
    """Test get_disease_genes aggregating genes across diseases."""

    @patch("drugdiscovery.databases.opentargets.post_with_retries")
    def test_collects_genes_from_diseases(self, mock_post):
        # Call 1: resolve target ID
        # Call 2: associated diseases (from get_target_diseases)
        # Call 3: disease targets for disease 1
        # Call 4: disease targets for disease 2
        mock_post.side_effect = [
            # _resolve_target_id call (inside get_target_diseases)
            _mock_graphql_response({
                "data": {
                    "search": {
                        "hits": [{"id": "ENSG00000146648"}]
                    }
                }
            }),
            # associated diseases query
            _mock_graphql_response({
                "data": {
                    "target": {
                        "id": "ENSG00000146648",
                        "approvedSymbol": "EGFR",
                        "associatedDiseases": {
                            "rows": [
                                {
                                    "disease": {
                                        "id": "EFO_0000311",
                                        "name": "lung carcinoma",
                                        "therapeuticAreas": [],
                                    },
                                    "score": 0.85,
                                },
                                {
                                    "disease": {
                                        "id": "EFO_0000616",
                                        "name": "glioblastoma",
                                        "therapeuticAreas": [],
                                    },
                                    "score": 0.72,
                                },
                            ]
                        },
                    }
                }
            }),
            # disease targets for lung carcinoma
            _mock_graphql_response({
                "data": {
                    "disease": {
                        "associatedTargets": {
                            "rows": [
                                {"target": {"approvedSymbol": "KRAS"}, "score": 0.9},
                                {"target": {"approvedSymbol": "TP53"}, "score": 0.8},
                                {"target": {"approvedSymbol": "EGFR"}, "score": 0.7},
                            ]
                        }
                    }
                }
            }),
            # disease targets for glioblastoma
            _mock_graphql_response({
                "data": {
                    "disease": {
                        "associatedTargets": {
                            "rows": [
                                {"target": {"approvedSymbol": "TP53"}, "score": 0.9},
                                {"target": {"approvedSymbol": "PTEN"}, "score": 0.85},
                                {"target": {"approvedSymbol": "IDH1"}, "score": 0.7},
                            ]
                        }
                    }
                }
            }),
        ]

        result = get_disease_genes("EGFR", limit=50)

        # Should contain genes from both diseases, deduplicated, without EGFR itself
        assert "KRAS" in result
        assert "TP53" in result
        assert "PTEN" in result
        assert "IDH1" in result
        assert "EGFR" not in result  # query gene excluded
        # Results are sorted
        assert result == sorted(result)

    @patch("drugdiscovery.databases.opentargets.post_with_retries")
    def test_returns_empty_when_no_diseases(self, mock_post):
        mock_post.side_effect = [
            _mock_graphql_response({
                "data": {
                    "search": {
                        "hits": [{"id": "ENSG00000000001"}]
                    }
                }
            }),
            _mock_graphql_response({
                "data": {
                    "target": {
                        "id": "ENSG00000000001",
                        "approvedSymbol": "OBSCURE",
                        "associatedDiseases": {
                            "rows": []
                        },
                    }
                }
            }),
        ]

        result = get_disease_genes("OBSCURE")
        assert result == []

    @patch("drugdiscovery.databases.opentargets.post_with_retries")
    def test_returns_empty_on_api_error(self, mock_post):
        mock_post.side_effect = Exception("Network error")

        result = get_disease_genes("EGFR")
        assert result == []

    @patch("drugdiscovery.databases.opentargets.post_with_retries")
    def test_handles_missing_disease_data(self, mock_post):
        mock_post.side_effect = [
            _mock_graphql_response({
                "data": {
                    "search": {
                        "hits": [{"id": "ENSG00000146648"}]
                    }
                }
            }),
            _mock_graphql_response({
                "data": {
                    "target": {
                        "id": "ENSG00000146648",
                        "approvedSymbol": "EGFR",
                        "associatedDiseases": {
                            "rows": [
                                {
                                    "disease": {
                                        "id": "EFO_0000311",
                                        "name": "lung carcinoma",
                                        "therapeuticAreas": [],
                                    },
                                    "score": 0.85,
                                },
                            ]
                        },
                    }
                }
            }),
            # disease query returns null disease node
            _mock_graphql_response({
                "data": {
                    "disease": None
                }
            }),
        ]

        result = get_disease_genes("EGFR")
        # Should gracefully return empty (disease node was None)
        assert result == []


# ---------------------------------------------------------------------------
# GraphQL Request Tests
# ---------------------------------------------------------------------------

class TestGraphqlRequest:
    """Test the internal _graphql_request helper."""

    @patch("drugdiscovery.databases.opentargets.post_with_retries")
    def test_sends_correct_payload(self, mock_post):
        from drugdiscovery.databases.opentargets import _graphql_request

        mock_post.return_value = _mock_graphql_response({"data": {}})

        _graphql_request("query Test { field }", {"key": "value"})

        mock_post.assert_called_once()
        call_kwargs = mock_post.call_args
        payload = call_kwargs.kwargs.get("json") or call_kwargs[1].get("json")
        assert payload["query"] == "query Test { field }"
        assert payload["variables"] == {"key": "value"}

    @patch("drugdiscovery.databases.opentargets.post_with_retries")
    def test_uses_correct_endpoint(self, mock_post):
        from drugdiscovery.databases.opentargets import (
            _graphql_request,
            OPENTARGETS_GRAPHQL,
        )

        mock_post.return_value = _mock_graphql_response({"data": {}})

        _graphql_request("query { field }", {})

        call_args = mock_post.call_args
        assert call_args[0][0] == OPENTARGETS_GRAPHQL

    @patch("drugdiscovery.databases.opentargets.post_with_retries")
    def test_propagates_exceptions(self, mock_post):
        from drugdiscovery.databases.opentargets import _graphql_request

        mock_post.side_effect = RuntimeError("POST failed after 3 attempts")

        try:
            _graphql_request("query { field }", {})
            assert False, "Expected RuntimeError to propagate"
        except RuntimeError as exc:
            assert "POST failed" in str(exc)
