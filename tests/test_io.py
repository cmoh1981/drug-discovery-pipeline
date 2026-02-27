"""Tests for drugdiscovery.utils.io module."""

import tempfile
from pathlib import Path

from drugdiscovery.utils.io import (
    read_csv, write_csv,
    read_fasta, write_fasta,
    read_json, write_json,
    extract_sequence_from_pdb,
    parse_pdb_binding_coords,
    compute_centroid,
)


class TestCSV:
    def test_roundtrip(self):
        rows = [
            {"id": "1", "name": "test", "value": "42"},
            {"id": "2", "name": "foo", "value": "99"},
        ]
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False, mode="w") as f:
            path = f.name
        write_csv(path, rows)
        loaded = read_csv(path)
        assert len(loaded) == 2
        assert loaded[0]["id"] == "1"
        assert loaded[1]["value"] == "99"

    def test_empty_rows(self):
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False, mode="w") as f:
            path = f.name
        write_csv(path, [])
        # Should not crash, just warn


class TestFASTA:
    def test_roundtrip(self):
        entries = [
            ("protein1 | test", "ACDEFGHIKLMNPQ"),
            ("protein2 | test2", "RSTVWY"),
        ]
        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, mode="w") as f:
            path = f.name
        write_fasta(path, entries)
        loaded = read_fasta(path)
        assert len(loaded) == 2
        assert loaded[0][0] == "protein1 | test"
        assert loaded[0][1] == "ACDEFGHIKLMNPQ"
        assert loaded[1][1] == "RSTVWY"

    def test_multiline_sequence(self):
        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, mode="w") as f:
            f.write(">test\n")
            f.write("ACDEF\n")
            f.write("GHIKL\n")
            path = f.name
        loaded = read_fasta(path)
        assert loaded[0][1] == "ACDEFGHIKL"


class TestJSON:
    def test_roundtrip(self):
        data = {"key": "value", "number": 42, "nested": {"a": [1, 2, 3]}}
        with tempfile.NamedTemporaryFile(suffix=".json", delete=False, mode="w") as f:
            path = f.name
        write_json(path, data)
        loaded = read_json(path)
        assert loaded["key"] == "value"
        assert loaded["number"] == 42
        assert loaded["nested"]["a"] == [1, 2, 3]


class TestPDBParsing:
    SAMPLE_PDB = """ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 50.00           N
ATOM      2  CA  ALA A   1       2.000   3.000   4.000  1.00 50.00           C
ATOM      3  N   GLY A   2       3.000   4.000   5.000  1.00 60.00           N
ATOM      4  CA  GLY A   2       4.000   5.000   6.000  1.00 60.00           C
ATOM      5  N   ASP A   3       5.000   6.000   7.000  1.00 70.00           N
ATOM      6  CA  ASP A   3       6.000   7.000   8.000  1.00 70.00           C
END"""

    def test_extract_sequence(self):
        seq = extract_sequence_from_pdb(self.SAMPLE_PDB, chain="A")
        assert seq == "AGD"

    def test_parse_binding_coords(self):
        coords = parse_pdb_binding_coords(self.SAMPLE_PDB, [1, 3], chain="A")
        assert len(coords) == 2
        assert coords[0]["resnum"] == 1
        assert coords[0]["x"] == 2.0  # CA atom
        assert coords[1]["resnum"] == 3

    def test_compute_centroid(self):
        coords = [
            {"x": 0.0, "y": 0.0, "z": 0.0},
            {"x": 10.0, "y": 10.0, "z": 10.0},
        ]
        cx, cy, cz = compute_centroid(coords)
        assert cx == 5.0
        assert cy == 5.0
        assert cz == 5.0

    def test_empty_centroid(self):
        assert compute_centroid([]) == (0.0, 0.0, 0.0)
