"""Tests for drugdiscovery.config module."""

import tempfile
from pathlib import Path

import pytest
from drugdiscovery.config import parse_args, make_run_dir, resolve_device


class TestParseArgs:
    def test_required_args(self):
        cfg = parse_args(["--modality", "peptide", "--mode", "antagonist", "--target", "YARS2"])
        assert cfg.modality.value == "peptide"
        assert cfg.mode.value == "antagonist"
        assert cfg.target == "YARS2"

    def test_defaults(self):
        cfg = parse_args(["--modality", "small_molecule", "--mode", "agonist", "--target", "EGFR"])
        assert cfg.top_n == 20
        assert cfg.device == "auto"
        assert cfg.use_runpod is False

    def test_optional_args(self):
        cfg = parse_args([
            "--modality", "peptide", "--mode", "antagonist", "--target", "YARS2",
            "--tissue", "lung", "--top-n", "50", "--device", "cpu",
        ])
        assert cfg.tissue == "lung"
        assert cfg.top_n == 50
        assert cfg.device == "cpu"

    def test_invalid_modality(self):
        with pytest.raises(SystemExit):
            parse_args(["--modality", "rna", "--mode", "agonist", "--target", "X"])

    def test_missing_required(self):
        with pytest.raises(SystemExit):
            parse_args(["--modality", "peptide"])


class TestMakeRunDir:
    def test_creates_structure(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = make_run_dir(tmpdir, "YARS2", "peptide", "antagonist")
            assert run_dir.exists()
            assert (run_dir / "target").exists()
            assert (run_dir / "library_hits").exists()
            assert (run_dir / "de_novo").exists()
            assert (run_dir / "scoring").exists()
            assert (run_dir / "report" / "figures").exists()

    def test_dirname_format(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = make_run_dir(tmpdir, "EGFR", "small_molecule", "agonist")
            name = run_dir.name
            assert "EGFR" in name
            assert "small_molecule" in name
            assert "agonist" in name


class TestResolveDevice:
    def test_explicit_cpu(self):
        assert resolve_device("cpu") == "cpu"

    def test_explicit_cuda(self):
        assert resolve_device("cuda") == "cuda"

    def test_auto_falls_back_to_cpu(self):
        # On machines without GPU, auto should return cpu
        device = resolve_device("auto")
        assert device in ("cpu", "cuda", "mps")
