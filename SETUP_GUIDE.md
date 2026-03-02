# Drug Discovery Pipeline - Setup Guide

## Overview

An automated computational drug discovery pipeline that screens databases, generates de novo candidates, predicts ADMET properties, scores and ranks compounds, and produces interactive HTML reports.

**Supports:** Small molecules and peptides
**Targets:** Any human protein (by gene name)
**Modes:** Agonist, antagonist, inhibitor, activator, modulator

## Quick Start

```bash
# 1. Install Python 3.10+ (tested on 3.14)
# 2. Install dependencies
pip install -e .

# 3. Run a pipeline
python -m drugdiscovery --target EGFR --modality small_molecule --mode antagonist
python -m drugdiscovery --target ESRRG --modality small_molecule --mode agonist
python -m drugdiscovery --target GLP1R --modality peptide --mode agonist
```

Results are saved to `results/<timestamp>_<target>_<modality>_<mode>/` with an interactive HTML dashboard at `report/final_report.html`.

## Requirements

### Core (Required)

```
Python >= 3.10
numpy >= 1.24.0
pandas >= 2.0.0
requests >= 2.31.0
biopython >= 1.83
rdkit >= 2023.9.1
pydantic >= 2.5.0
jinja2 >= 3.1.0
tqdm >= 4.66.0
matplotlib >= 3.8.0
seaborn >= 0.13.0
pyyaml >= 6.0
python-dotenv >= 1.0.0
```

### Optional

| Extra | Install | Purpose |
|-------|---------|---------|
| GPU | `pip install -e ".[gpu]"` | PepMLM peptide generation (torch, transformers) |
| Docking | `pip install -e ".[docking]"` | AutoDock Vina molecular docking |
| Cloud | `pip install -e ".[cloud]"` | RunPod GPU dispatch |
| Web | `pip install -e ".[web]"` | FastAPI SaaS service layer |
| Dev | `pip install -e ".[dev]"` | pytest, coverage, httpx |
| All | `pip install -e ".[all]"` | Everything above |

### RDKit Installation

RDKit is the most critical dependency. Install via:

```bash
# Option 1: pip (recommended)
pip install rdkit

# Option 2: conda
conda install -c conda-forge rdkit
```

## Installation Steps

### 1. Clone or Unzip

```bash
# From zip file
unzip drug-discovery-pipeline.zip
cd drug-discovery-pipeline

# Or from git
git clone <repo-url>
cd drug-discovery-pipeline
```

### 2. Create Virtual Environment (Recommended)

```bash
python -m venv .venv

# Windows
.venv\Scripts\activate

# macOS/Linux
source .venv/bin/activate
```

### 3. Install

```bash
# Core only (runs pipeline without GPU/docking)
pip install -e .

# With dev tools
pip install -e ".[dev]"

# Everything
pip install -e ".[all]"
```

### 4. Verify Installation

```bash
# Run tests (should show 224 passed)
pytest tests/ -v

# Quick pipeline test
python -m drugdiscovery --target EGFR --modality small_molecule --mode antagonist
```

## Usage

### Command Line

```bash
python -m drugdiscovery \
    --target <GENE_NAME> \
    --modality <small_molecule|peptide> \
    --mode <agonist|antagonist|inhibitor|activator|modulator>
```

### Examples

```bash
# Small molecule antagonist for EGFR (cancer target)
python -m drugdiscovery --target EGFR --modality small_molecule --mode antagonist

# Peptide agonist for GLP1R (diabetes target)
python -m drugdiscovery --target GLP1R --modality peptide --mode agonist

# Small molecule agonist for ESRRG (metabolic target)
python -m drugdiscovery --target ESRRG --modality small_molecule --mode agonist
```

### Custom Configuration

```bash
python -m drugdiscovery --target BRAF --modality small_molecule --mode antagonist \
    --config configs/defaults.yaml
```

Edit `configs/defaults.yaml` to customize scoring weights, database sources, docking parameters, etc.

## Pipeline Modules

The pipeline runs 10 modules sequentially:

| Module | Name | Description |
|--------|------|-------------|
| M1 | Target Preparation | Resolves gene to UniProt, fetches sequence, identifies anti-targets |
| M2 | Library Screening | Queries ChEMBL, ZINC22/PubChem, BindingDB, DrugBank, OpenTargets, COCONUT |
| M3 | De Novo Design | BRICS retrosynthetic recombination from 24 FDA-approved drug scaffolds |
| M4 | Structure Prediction | AlphaFold structure download or ColabFold/Boltz prediction |
| M4.5 | Molecular Docking | AutoDock Vina / Gnina / DiffDock (when PDB available) |
| M4.6 | Perturbation Biology | CMAP connectivity scoring + STRING PPI network effects |
| M7 | ADMET Prediction | ADMETlab3 API or RDKit fallback (solubility, toxicity, metabolism) |
| M5 | Scoring & Ranking | Composite scoring: binding + selectivity + drug-likeness + ADMET + perturbation |
| M8 | Delivery System | Route-of-administration recommendations based on tissue target |
| M9 | Final Report | Interactive HTML dashboard with Chart.js + CSV exports |

## Output Structure

```
results/<timestamp>_<target>_<modality>_<mode>/
  config_snapshot.yaml          # Run configuration
  pipeline.log                  # Full execution log
  target/
    target_profile.json         # UniProt data, sequence, anti-targets
    target.fasta                # Protein sequence
    binding_pockets.json        # Detected binding pockets
    anti_targets.json           # Related proteins to avoid
  library_hits/
    chembl_hits.csv             # ChEMBL database hits
    zinc22_hits.csv             # ZINC22/PubChem hits
    merged_hits.csv             # Deduplicated merged hits
  de_novo/
    candidates.csv              # BRICS-generated novel molecules
  structure/                    # AlphaFold PDB files
  docking/                      # Docking poses and scores
  perturbation/
    disease_genes.json          # Disease-associated genes
    network_effect.json         # PPI network analysis
  admet/
    admet_profiles.csv          # Full ADMET predictions
  scoring/
    scored_candidates.csv       # All scored candidates
    top_20_candidates.csv       # Top 20 by composite score
  delivery/
    delivery_recommendations.csv
  report/
    final_report.html           # Interactive HTML dashboard
    final_report.md             # Markdown report
    executive_summary.csv       # Top candidates summary
    figures/                    # Distribution plots (PNG)
```

## Database Sources

| Database | Data Type | Access | Status |
|----------|-----------|--------|--------|
| ChEMBL | Bioactivity data | Free API | Working |
| ZINC22/PubChem | Purchasable compounds | Free API | Working |
| OpenTargets | Known drugs | Free API | Working |
| BindingDB | Binding affinities | Free API | SSL expired (circuit breaker skips) |
| DrugBank | Drug data | Requires account | Needs download |
| COCONUT | Natural products | Free API | API endpoint changed |
| ADMETlab3 | ADMET predictions | Free API | SSL expired (falls back to RDKit) |

### Setting Up DrugBank (Optional)

1. Create a free academic account at https://go.drugbank.com
2. Download the full database:
   ```bash
   curl -Lfv -o drugbank_full.zip -u YOUR_EMAIL:YOUR_PASSWORD \
       https://go.drugbank.com/releases/5-1-14/downloads/all-full-database
   mkdir -p data/drugbank
   unzip drugbank_full.zip -d data/drugbank/
   ```

### Setting Up Local BindingDB (Optional)

```bash
git clone --depth 1 https://github.com/dhimmel/bindingdb.git data/bindingdb
```

## Scoring Weights

Default composite score weights (configurable in `configs/defaults.yaml`):

| Component | Weight | Description |
|-----------|--------|-------------|
| Binding Energy | 0.30 | Docking score or predicted affinity |
| Selectivity | 0.20 | Target vs anti-target selectivity |
| ADMET | 0.15 | Aggregate ADMET safety score |
| Perturbation | 0.15 | CMAP connectivity + PPI network effects |
| Drug-likeness | 0.10 | Lipinski/Veber rule compliance |
| MoA Consistency | 0.10 | Mode-of-action prediction match |

## Troubleshooting

### "No module named rdkit"
```bash
pip install rdkit
# or: conda install -c conda-forge rdkit
```

### SSL Certificate Errors
External APIs (BindingDB, ADMETlab3) occasionally have expired SSL certs. The pipeline has circuit breakers that automatically skip after the first failure and fall back to local methods.

### Slow Pipeline Run
If M5 scoring is slow, it's likely hitting an external API with connection issues. The circuit breaker should trip after the first failure. If running an older version, update from this zip.

### Unicode Errors on Windows
The pipeline handles Windows cp949/cp1252 encoding. If you see encoding errors in logs, ensure your terminal supports UTF-8:
```bash
chcp 65001
```

### Tests Failing
```bash
# Run full test suite
pytest tests/ -v

# Expected: 224 passed
```

## Project Structure

```
drug-discovery-pipeline/
  drugdiscovery/              # Core pipeline package
    __init__.py
    __main__.py               # CLI entry point
    config.py                 # Configuration management
    pipeline.py               # Main pipeline orchestrator
    types.py                  # Data types (Candidate, PipelineConfig, etc.)
    databases/                # Database connectors
      chembl.py               # ChEMBL bioactivity database
      zinc22.py               # ZINC22 via PubChem similarity
      bindingdb.py            # BindingDB binding affinities
      drugbank.py             # DrugBank drug database
      opentargets.py          # Open Targets known drugs
      coconut.py              # COCONUT natural products
      pubchem_db.py           # PubChem bioassay data
      peptide_dbs.py          # Peptide-specific databases
    modules/                  # Pipeline modules (M1-M9)
      target_prep.py          # M1: Target preparation
      library_screening.py    # M2: Database screening
      de_novo_sm.py           # M3: Small molecule de novo design (BRICS)
      de_novo_peptide.py      # M3: Peptide de novo design (PepMLM)
      structure_pred.py       # M4: Structure prediction
      docking.py              # M4.5: Molecular docking
      perturbation.py         # M4.6: Perturbation biology
      admet.py                # M7: ADMET prediction
      scoring.py              # M5: Composite scoring
      delivery.py             # M8: Delivery system
      reporting.py            # M9: HTML report generation
      modifications.py        # Chemical modifications
    tools/                    # External tool integrations
      admetlab3.py            # ADMETlab3 REST API
      vina.py                 # AutoDock Vina
      gnina.py                # Gnina CNN docking
      diffdock.py             # DiffDock ML docking
      boltz.py                # Boltz structure prediction
      colabfold.py            # ColabFold integration
      pepmlm.py               # PepMLM peptide generation
      pubchem.py              # PubChem utilities
      cmap_client.py          # CMAP connectivity client
      string_ppi.py           # STRING PPI network
      pocket_detection.py     # Binding pocket detection
      deepdta.py              # DeepDTA binding affinity
      posecheck_tool.py       # Pose quality checking
      diffpepbuilder.py       # DiffPepBuilder integration
    utils/                    # Shared utilities
      web.py                  # HTTP client with retries
      io.py                   # File I/O helpers
      chemistry.py            # Chemistry utilities
      compute.py              # RunPod GPU dispatch
  service/                    # FastAPI web service (optional)
  tests/                      # Unit tests (224 tests)
  tests_service/              # Service layer tests
  configs/
    defaults.yaml             # Default configuration
  pyproject.toml              # Package configuration
  requirements.txt            # Pip requirements
```

## License

MIT
