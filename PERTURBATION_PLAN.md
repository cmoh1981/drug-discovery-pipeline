# Perturbation Biology Integration Plan

## Executive Summary

Integrate perturbation biology scoring into the drug discovery pipeline as module **M4.6**
(between Docking and ADMET). This adds a **Perturbation Response Score** that predicts
how each drug candidate will perturb disease-relevant gene expression at the cellular level.

### Why This Matters

Traditional drug discovery scores candidates on binding affinity and physicochemical
properties alone. Perturbation biology adds a critical dimension: **will this compound
actually produce the desired cellular effect?** A drug can bind perfectly but fail to
produce the intended transcriptional response — perturbation scoring catches this.

---

## Research Synthesis

### Key Findings from 60+ Models Reviewed

**1. Nature Methods 2025 (Ahlmann-Eltze et al.)**
- Deep learning models (scGPT, Geneformer, scFoundation, etc.) do NOT outperform
  simple linear baselines for gene perturbation prediction
- Simple mean-shift baseline is surprisingly competitive
- Implication: **Don't overbuild. Simple methods with good embeddings win.**

**2. GenBio AI / Nature Methods 2025 Benchmark (Wei et al.)**
- Interactome-based embeddings (WaveGC, GenePT using STRING PPI) consistently
  outperform sequence-based models
- kNN on good embeddings beats diffusion, flow matching, and transformers
- Gene-gene interaction networks are more informative than gene sequences

**3. LPM — Large Perturbation Model (Nature Computational Science 2025)**
- Integrates BOTH chemical and genetic perturbations in a unified framework
- Disentangles perturbation, readout, and context as separate dimensions
- Identifies shared mechanisms of action between drugs and gene knockouts
- Best existing model for cross-context perturbation prediction

**4. PRnet (Nature Communications 2024)**
- Deep generative model for chemical perturbation response
- Predicts transcriptional responses to novel (unseen) compounds
- Uses molecular structure + cell context as input
- 72 stars, 10 forks on GitHub — active community

**5. PerturbNet (Molecular Systems Biology 2025)**
- Predicts single-cell responses to unseen chemical AND genetic perturbations
- Accounts for dosage and cell type covariates
- Chemical structure → gene expression change prediction

**6. CPA — Compositional Perturbation Autoencoder (Theis Lab)**
- OOD prediction of unseen drug combinations
- Interpretable drug embeddings
- Dose-response curve estimation
- Well-maintained codebase (99 stars, BSD-3)

**7. GEARS (Nature Biotechnology 2023)**
- Graph-enhanced gene activation/repression simulator
- Uses knowledge graph of gene-gene relationships
- 40% higher precision for genetic interaction prediction
- Primarily for genetic (not chemical) perturbations

**8. Connectivity Map / L1000 (Broad Institute)**
- 1.3M gene expression profiles across ~27K compounds
- Pre-computed drug signatures — no ML needed
- Connectivity scoring: does drug reverse disease signature?
- Most practical, battle-tested approach for drug discovery

**9. scDrugMap Benchmark (Nature Communications 2025)**
- Benchmarks 8 foundation models + 2 LLMs across 495K cells
- scFoundation best in pooled evaluation
- UCE best after fine-tuning
- scGPT best in zero-shot settings

**10. Systema Framework (Nature Biotechnology 2025)**
- Shows current metrics overestimate perturbation prediction performance
- "Systematic variation" confounds evaluation
- Most methods learn selection bias, not true perturbation effects

---

## Architecture Decision

### Chosen Approach: Hybrid L1000/CMAP + Lightweight Embedding Model

Based on the research, we will NOT build a complex deep learning perturbation predictor.
Instead, we use a **pragmatic three-tier approach**:

#### Tier 1: L1000/CMAP Connectivity Score (Primary — No GPU Needed)
- Query Broad's CLUE API or use pre-downloaded L1000 signatures
- For each candidate SMILES/name, find matching or similar compound signatures
- Compute connectivity score: does the candidate's signature reverse the disease gene set?
- **This is the most validated approach in real drug discovery**

#### Tier 2: Perturbation Embedding Similarity (Secondary — CPU)
- Use pre-computed drug perturbation embeddings from ChemicalProbes/LINCS
- For novel compounds: compute Morgan fingerprint → kNN to nearest L1000 compound
- Transfer the perturbation profile from the nearest neighbor
- Based on the finding that kNN on good embeddings beats complex models

#### Tier 3: Target Gene Network Effect (Tertiary — CPU)
- Use STRING PPI network to model downstream effects of target inhibition/activation
- For each candidate's target: propagate effect through interactome
- Score: what fraction of disease-relevant genes are affected?
- Based on finding that interactome embeddings outperform all others

### What We Will NOT Do
- Train our own foundation model (waste of compute, beaten by baselines)
- Fine-tune scGPT/Geneformer (unreliable per Nature Methods 2025)
- Build complex VAE/diffusion models (overkill for scoring, CPA-level complexity)
- Require GPU for perturbation scoring (must work on CPU for accessibility)

---

## Implementation Plan

### New Files to Create

```
drugdiscovery/
├── modules/
│   └── perturbation.py          # M4.6 module — main entry point
├── tools/
│   ├── cmap_client.py           # L1000/CMAP API client
│   └── string_ppi.py            # STRING PPI network effects
```

### Files to Modify

```
drugdiscovery/
├── types.py                     # Add perturbation fields to Candidate
├── pipeline.py                  # Wire M4.6 between M4.5 and M7
├── modules/scoring.py           # Add perturbation_score to composite
├── modules/reporting.py         # Include perturbation data in reports
├── configs/defaults.yaml        # Add perturbation config + scoring weight
```

### Step 1: Extend Candidate dataclass (types.py)

Add these fields to `Candidate`:
```python
perturbation_score: float = 0.0       # 0-1 aggregate perturbation score
cmap_connectivity: float = 0.0        # L1000/CMAP connectivity score
cmap_compound_match: str = ""         # Nearest CMAP compound matched
network_effect_score: float = 0.0     # STRING PPI network propagation score
disease_signature_reversal: float = 0.0  # How well drug reverses disease sig
perturbation_metadata: dict = field(default_factory=dict)
```

Update `to_dict()` to include perturbation fields.

### Step 2: Build L1000/CMAP Client (tools/cmap_client.py)

```python
def query_cmap_signature(smiles: str, gene_name: str = "") -> dict | None:
    """Query CLUE API for compound perturbation signature.

    Strategy:
    1. Try exact SMILES match in LINCS compounds
    2. Try gene name if compound is a known inhibitor
    3. Fall back to Tanimoto-nearest neighbor in L1000 library

    Returns dict with:
      - connectivity_score: float (-1 to 1)
      - matched_compound: str
      - up_genes: list[str]  (top 50 upregulated)
      - down_genes: list[str]  (top 50 downregulated)
      - cell_lines: list[str]
    """
```

### Step 3: Build STRING PPI Network Tool (tools/string_ppi.py)

```python
def compute_network_effect(
    target_gene: str,
    disease_genes: list[str],
    mode: str = "antagonist",
    depth: int = 2,
) -> dict:
    """Compute network propagation effect score via STRING API.

    1. Get first-order interactors of target_gene from STRING
    2. Propagate perturbation effect through network (depth hops)
    3. Score overlap with disease-relevant genes

    Returns dict with:
      - network_effect_score: float (0-1)
      - affected_disease_genes: list[str]
      - pathway_enrichment: list[dict]
    """
```

### Step 4: Build Perturbation Module (modules/perturbation.py)

```python
def predict_perturbation(
    cfg: PipelineConfig,
    candidates: list[Candidate],
    target_profile: TargetProfile,
    output_dir: Path,
) -> list[Candidate]:
    """M4.6 — Predict cellular perturbation effects for all candidates.

    Three-tier scoring:
      T1: L1000/CMAP connectivity (weight 0.50)
      T2: Tanimoto-kNN perturbation transfer (weight 0.20)
      T3: STRING PPI network effect (weight 0.30)

    Sets candidate.perturbation_score (0-1, higher = better reversal of
    disease signature / stronger on-target network effect).
    """
```

### Step 5: Wire into Pipeline (pipeline.py)

Insert M4.6 between M4.5 (Docking) and M7 (ADMET):

```python
# --- M4.6: Perturbation Biology ---
logger.info("[M4.6] Perturbation Biology Scoring")
try:
    from drugdiscovery.modules.perturbation import predict_perturbation
    candidates = predict_perturbation(cfg, candidates, target_profile, run_dir / "perturbation")
    logger.info("[M4.6] Perturbation scores computed for %d candidates", len(candidates))
except Exception as exc:
    logger.warning("[M4.6] Perturbation scoring skipped: %s", exc, exc_info=True)
```

### Step 6: Update Composite Scoring (scoring.py)

Add `perturbation` to the composite score components:

```python
# New weights (rebalanced):
scoring_weights = {
    "binding_energy": 0.30,     # was 0.35
    "selectivity": 0.20,        # was 0.25
    "drug_likeness": 0.10,      # was 0.15
    "admet_aggregate": 0.15,    # unchanged
    "moa_consistency": 0.10,    # unchanged
    "perturbation": 0.15,       # NEW
}
```

### Step 7: Update Config (defaults.yaml)

```yaml
perturbation:
  enable: true
  cmap_api_url: "https://api.clue.io/api"
  cmap_api_key: ""               # User provides their CLUE API key
  string_api_url: "https://string-db.org/api"
  disease_genes: []              # Auto-populated from target_prep
  network_depth: 2
  knn_neighbors: 5
  tier_weights:
    cmap_connectivity: 0.50
    knn_transfer: 0.20
    network_effect: 0.30
```

### Step 8: Update Reporting (reporting.py)

Add perturbation section to final report:
- Perturbation score distribution
- Top CMAP compound matches
- Network effect visualization (genes affected)
- Disease signature reversal assessment

### Step 9: Tests

Add tests in `tests/test_perturbation.py`:
- Unit tests for each tier of scoring
- Integration test for full module
- Mock API responses for CMAP and STRING
- Edge cases (no SMILES, no matches, API failures)

---

## Disease Gene Set Derivation

The perturbation module needs disease-relevant genes to score against.
These are derived automatically from the target:

1. **From target_prep (M1)**: GO biological process terms → gene sets
2. **From STRING PPI**: First-order interactors of target gene
3. **From KEGG pathways**: Pathways containing the target gene
4. **Hardcoded disease gene sets**: For common targets (EGFR, BRAF, etc.)

---

## Graceful Degradation

Every tier has fallbacks:

| Tier | Primary | Fallback | Default |
|------|---------|----------|---------|
| T1 CMAP | CLUE API query | Pre-cached L1000 data | 0.5 |
| T2 kNN | Morgan FP similarity to L1000 | SMILES substring match | 0.5 |
| T3 Network | STRING API query | Cached PPI for common targets | 0.5 |

If ALL tiers fail → `perturbation_score = 0.5` (neutral, doesn't help or hurt).

---

## Expected Impact

- Candidates that bind well BUT produce wrong cellular effects get penalized
- Candidates matching known drug signatures for the disease get boosted
- Network-connected candidates (hitting right pathways) score higher
- More biologically meaningful rankings than binding + ADMET alone
