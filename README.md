# Drug Discovery Pipeline

Unified drug discovery platform supporting both small molecules and peptides, with agonist/antagonist targeting, library screening, de novo design, ADMET prediction, and delivery system recommendations.

## Quick Start

```bash
# Install
pip install -e .

# Run peptide antagonist pipeline
python -m drugdiscovery --modality peptide --mode antagonist --target YARS2

# Run small molecule agonist pipeline
python -m drugdiscovery --modality small_molecule --mode agonist --target EGFR
```

## Features

- **Dual modality**: Small molecules and peptides in one unified pipeline
- **9-module architecture**: Target prep → Library screening → De novo design → Structure prediction → Scoring → Modifications → ADMET → Delivery → Reporting
- **Multi-database screening**: ChEMBL, DrugBank, COCONUT, ZINC22, BIOPEP, APD3, DBAASP
- **De novo design**: PepMLM + DiffPepBuilder (peptides), REINVENT4/RDKit (small molecules)
- **AI scoring**: Composite scoring with binding, selectivity, drug-likeness, diversity
- **ADMET prediction**: Peptide-specific local computation + ADMETlab 3.0 API for small molecules
- **Delivery recommendations**: Tissue-specific delivery system mapping
- **GPU support**: Auto-detection with CPU fallback, RunPod cloud GPU integration

## CLI Arguments

| Argument | Required | Values | Description |
|----------|----------|--------|-------------|
| `--modality` | Yes | `small_molecule`, `peptide` | Drug modality |
| `--mode` | Yes | `agonist`, `antagonist` | Mode of action |
| `--target` | Yes | Gene name, UniProt ID, or PDB ID | Target protein |
| `--tissue` | No | e.g., `lung`, `liver`, `brain` | Target tissue (auto-detected if omitted) |
| `--top-n` | No | Integer (default: 20) | Top candidates to carry forward |
| `--device` | No | `auto`, `cpu`, `cuda` | Compute device |
| `--use-runpod` | No | Flag | Use RunPod for GPU tasks |
| `--output` | No | Path (default: `results/`) | Output directory |
| `--config` | No | YAML path | Config overrides |

## Pipeline Modules

| Module | Function | Tools |
|--------|----------|-------|
| M1: Target Prep | Resolve target, fetch structure, detect pockets | UniProt API, AlphaFold DB |
| M2: Library Screening | Multi-database search and merge | ChEMBL, DrugBank, COCONUT, ZINC22 |
| M3: De Novo Design | Generate novel candidates | PepMLM, DiffPepBuilder, RDKit |
| M4: Structure Prediction | 3D structure generation | ESMFold, ColabFold, RDKit ETKDG |
| M5: Scoring | Composite binding & property scoring | Physics-based + AI |
| M6: Modifications | Chemical modifications | Cyclization, stapling, PEGylation, etc. |
| M7: ADMET | Absorption, distribution, metabolism, excretion, toxicity | Local + ADMETlab 3.0 |
| M8: Delivery | Tissue-specific delivery recommendation | Rule-based engine |
| M9: Reporting | HTML + Markdown + CSV reports | Matplotlib, Jinja2 |

## Output Structure

Each run creates a timestamped directory:

```
results/<timestamp>_<target>_<modality>_<mode>/
    config_snapshot.yaml
    pipeline.log
    target/                  # Target profile, FASTA, PDB, pockets
    library_hits/            # Per-database CSV hits + merged
    de_novo/                 # Generated candidates
    scoring/                 # Scored and ranked candidates
    modifications/           # Modified variants
    admet/                   # ADMET profiles
    delivery/                # Delivery recommendations
    report/                  # HTML, Markdown, CSV, figures
```

## Configuration

Override defaults via YAML config:

```yaml
scoring_weights:
  binding_energy: 0.30
  structure_confidence: 0.20
  selectivity: 0.20
  drug_likeness: 0.15
  sequence_diversity: 0.15

pepmlm:
  lengths: [8, 10, 12, 15, 20]
  num_per_length: 100
  top_k: 5
```

## Requirements

Core: `numpy`, `pandas`, `requests`, `biopython`, `rdkit`, `pydantic`, `jinja2`, `tqdm`, `matplotlib`, `seaborn`, `pyyaml`

Optional: `torch`, `transformers` (GPU), `vina`, `meeko` (docking), `runpod` (cloud GPU)

```bash
pip install -e ".[all]"   # Install everything
pip install -e ".[gpu]"   # GPU support only
```

## Testing

```bash
pip install -e ".[dev]"
pytest tests/ -v
```

## Project Structure

```
drugdiscovery/
    __main__.py          # CLI entry point
    pipeline.py          # Orchestrator
    config.py            # Config + CLI parsing
    types.py             # Data types
    modules/             # 9 pipeline modules (M1-M9)
    databases/           # Database API clients
    tools/               # External tool wrappers
    utils/               # Shared utilities
```
