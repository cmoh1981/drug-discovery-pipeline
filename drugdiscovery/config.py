"""Configuration parsing and CLI argument handling."""

from __future__ import annotations

import argparse
import logging
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

import yaml

from drugdiscovery.types import Modality, ModeOfAction, PipelineConfig

logger = logging.getLogger(__name__)


def parse_args(argv: list[str] | None = None) -> PipelineConfig:
    """Parse CLI arguments and return a PipelineConfig."""
    parser = argparse.ArgumentParser(
        prog="drugdiscovery",
        description="Unified drug discovery pipeline for small molecules and peptides",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--modality",
        required=True,
        choices=["small_molecule", "peptide"],
        help="Drug modality",
    )
    parser.add_argument(
        "--mode",
        required=True,
        choices=["agonist", "antagonist"],
        help="Mode of action",
    )
    parser.add_argument(
        "--target",
        required=True,
        help="Target identifier: gene name, UniProt ID, or PDB ID",
    )
    parser.add_argument(
        "--tissue",
        default="",
        help="Target tissue for delivery (auto-detected if omitted)",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=20,
        help="Number of top candidates to carry forward",
    )
    parser.add_argument(
        "--device",
        default="auto",
        choices=["auto", "cpu", "cuda"],
        help="Compute device",
    )
    parser.add_argument(
        "--use-runpod",
        action="store_true",
        help="Use RunPod for GPU tasks",
    )
    parser.add_argument(
        "--output",
        default="results/",
        help="Output directory",
    )
    parser.add_argument(
        "--config",
        default="",
        help="Path to YAML config file for advanced overrides",
    )

    args = parser.parse_args(argv)

    cfg = PipelineConfig(
        modality=Modality(args.modality),
        mode=ModeOfAction(args.mode),
        target=args.target,
        tissue=args.tissue,
        top_n=args.top_n,
        device=args.device,
        use_runpod=args.use_runpod,
        output_dir=args.output,
        config_file=args.config,
    )

    if args.config and Path(args.config).is_file():
        cfg = _apply_yaml_overrides(cfg, args.config)

    return cfg


def _apply_yaml_overrides(cfg: PipelineConfig, yaml_path: str) -> PipelineConfig:
    """Apply overrides from a YAML config file."""
    with open(yaml_path, "r") as fh:
        overrides: dict[str, Any] = yaml.safe_load(fh) or {}

    if "target_pdb_path" in overrides:
        cfg.target_pdb_path = str(overrides["target_pdb_path"])
    if "scoring_weights" in overrides:
        sw = overrides["scoring_weights"]
        if not isinstance(sw, dict):
            raise TypeError(f"scoring_weights must be a dict, got {type(sw).__name__}")
        for k, v in sw.items():
            if not isinstance(v, (int, float)):
                raise TypeError(f"scoring_weights['{k}'] must be numeric, got {type(v).__name__}")
        cfg.scoring_weights.update(sw)
    if "pepmlm" in overrides:
        if not isinstance(overrides["pepmlm"], dict):
            raise TypeError(f"pepmlm must be a dict, got {type(overrides['pepmlm']).__name__}")
        cfg.pepmlm_settings.update(overrides["pepmlm"])
    if "diffpepbuilder" in overrides:
        if not isinstance(overrides["diffpepbuilder"], dict):
            raise TypeError(f"diffpepbuilder must be a dict, got {type(overrides['diffpepbuilder']).__name__}")
        cfg.diffpepbuilder_settings.update(overrides["diffpepbuilder"])
    if "vina" in overrides:
        if not isinstance(overrides["vina"], dict):
            raise TypeError(f"vina must be a dict, got {type(overrides['vina']).__name__}")
        cfg.vina_settings.update(overrides["vina"])
    if "library" in overrides:
        if not isinstance(overrides["library"], dict):
            raise TypeError(f"library must be a dict, got {type(overrides['library']).__name__}")
        cfg.library_settings.update(overrides["library"])

    return cfg


def resolve_device(device: str) -> str:
    """Resolve 'auto' device to actual device string."""
    if device != "auto":
        return device
    try:
        import torch
        if torch.cuda.is_available():
            return "cuda"
    except ImportError:
        pass
    return "cpu"


def make_run_dir(output_dir: str, target: str, modality: str, mode: str) -> Path:
    """Create a timestamped run output directory."""
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    name = f"{ts}_{target}_{modality}_{mode}"
    run_dir = Path(output_dir) / name
    subdirs = [
        "target", "library_hits", "de_novo",
        "scoring", "modifications", "admet",
        "delivery", "report/figures",
    ]
    for sub in subdirs:
        (run_dir / sub).mkdir(parents=True, exist_ok=True)
    return run_dir


def setup_logging(run_dir: Path) -> None:
    """Configure logging to both console and file."""
    log_file = run_dir / "pipeline.log"
    root = logging.getLogger()
    root.setLevel(logging.INFO)

    # Prevent duplicate handlers on repeated calls
    root.handlers.clear()

    fmt = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    fh = logging.FileHandler(str(log_file), encoding="utf-8")
    fh.setLevel(logging.INFO)
    fh.setFormatter(fmt)
    root.addHandler(fh)

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    root.addHandler(ch)


def save_config_snapshot(cfg: PipelineConfig, run_dir: Path) -> None:
    """Save frozen config to run directory."""
    snapshot = {
        "modality": cfg.modality.value,
        "mode": cfg.mode.value,
        "target": cfg.target,
        "tissue": cfg.tissue,
        "top_n": cfg.top_n,
        "device": cfg.device,
        "use_runpod": cfg.use_runpod,
        "scoring_weights": cfg.scoring_weights,
        "pepmlm": cfg.pepmlm_settings,
        "diffpepbuilder": cfg.diffpepbuilder_settings,
        "vina": cfg.vina_settings,
        "library": cfg.library_settings,
    }
    path = run_dir / "config_snapshot.yaml"
    with open(path, "w") as fh:
        yaml.dump(snapshot, fh, default_flow_style=False, sort_keys=False)
    logger.info("Config snapshot saved to %s", path)
