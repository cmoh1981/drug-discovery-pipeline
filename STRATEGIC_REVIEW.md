# Drug Discovery Pipeline — Strategic Review & Improvement Plan

**Date:** 2026-02-27
**Scope:** Full architectural review + state-of-the-art benchmarking

---

## Executive Summary

The pipeline has **excellent engineering foundations** (92 tests, graceful degradation, dual modality) but has **critical scientific gaps** that prevent it from producing commercially viable candidates. The two most urgent issues are:

1. **ADMET is not included in composite scoring** — candidates with fatal pharmacological liabilities rank alongside clean ones
2. **Selectivity scoring is scientifically unsound** — amino acid composition similarity ≠ binding selectivity

The competitive landscape has shifted dramatically. Leading platforms (Insilico, Recursion, Iambic) now operate **closed-loop design-make-test-learn cycles in under one week**, with programs reaching Phase 2-3 trials. Our pipeline needs both scientific and business upgrades to be competitive.

---

## Part 1: Critical Scientific Defects (Fix Immediately)

### 1.1 ADMET Not in Composite Score

**Current:** `composite_score = binding(0.30) + structure_conf(0.20) + selectivity(0.20) + drug_likeness(0.15) + diversity(0.15)`

ADMET (M7) runs AFTER ranking decisions are made. A candidate with hepatotoxicity, zero solubility, and 5 ADMET flags ranks identically to a clean candidate.

**Fix:** Add `admet_aggregate` to composite scoring with weight 0.15. Reorder pipeline: M7 (ADMET) before M5 (final scoring).

### 1.2 Scoring Methodology Issues

| Problem | Current Behavior | Fix |
|---------|-----------------|-----|
| `sequence_diversity` in score | Portfolio metric treated as quality metric; same peptide gets different scores in different pools | Remove from composite; apply as post-ranking filter |
| `structure_confidence` as weight | High pLDDT ≠ good drug; rewards well-predicted poor binders | Use as gate (exclude pLDDT < 0.5) not weighted component |
| Batch-relative binding normalization | Best binder always gets 1.0 even if -3.0 kcal/mol (weak) | Use absolute thresholds (discard > -5.0 kcal/mol) |
| Modified candidates not re-docked | PEGylated peptide (+2000-5000 Da) keeps parent binding score | Re-dock after M6 modifications |

**Proposed revised scoring:**
```yaml
scoring_weights:
  binding_energy:     0.35  # absolute-thresholded
  selectivity:        0.25  # docking-based, not composition
  drug_likeness:      0.15
  admet_aggregate:    0.15
  moa_consistency:    0.10
# structure_confidence: gate at > 0.5 (not scored)
# sequence_diversity: post-ranking filter (not scored)
```

### 1.3 Dead Code / Stubs Masquerading as Features

| Component | Status | Impact |
|-----------|--------|--------|
| DeepDTA (`deepdta.py`) | Always returns `None`; never imported by any module | No binding prediction beyond Vina |
| ADMETlab 3.0 client (`admetlab3.py`) | Complete API client, never called from `admet.py` | SM ADMET uses only basic RDKit rules |
| REINVENT 4 (`de_novo_sm.py`) | Stub that always returns `[]` | SM de novo = random fragment assembly |
| SM selectivity | Flat `0.5` placeholder | No selectivity computation for small molecules |

---

## Part 2: State-of-the-Art Tool Upgrades

### 2.1 Peptide Design (M3) — Replace/Supplement PepMLM + DiffPepBuilder

| Tool | What It Does | Why Better | Status |
|------|-------------|-----------|--------|
| **RFpeptides** (Baker Lab, 2024) | Macrocyclic peptide design via RFdiffusion | Validated 6 nM KD (GABARAP); handles ring-closure natively; X-ray confirmed | Open source (RoseTTAFold) |
| **PepCCD** | Target-conditioned peptide diffusion (sequence-only input) | Works without crystal structure; outperforms sequence baselines | Available |
| **ProteinMPNN** | Sequence design given backbone structure | Standard for fixed-backbone redesign | Open source |
| **AfCycDesign** | AlphaFold2 with cyclic offset for cyclic peptide prediction | Specifically designed for cyclic peptides | Open source |

### 2.2 Structure Prediction (M4) — Beyond ColabFold/ESMFold

| Tool | Advantage Over Current | Access |
|------|----------------------|--------|
| **AlphaFold3** (EBI server) | Predicts protein-ligand complexes directly; 70-80% success for peptide-protein | Free web server |
| **Boltz-1** | Open-source AF3 equivalent; fully local | MIT license, GitHub |
| **Boltz-2** | Structure + binding affinity jointly predicted; 1000x faster than FEP at comparable accuracy | MIT license |
| **Chai-1** | Multi-modal: proteins, small molecules, DNA, RNA, glycosylations | Partially open |

### 2.3 Docking & Scoring (M4.5/M5) — Beyond AutoDock Vina

| Tool | Improvement | Access |
|------|------------|--------|
| **GNINA** | CNN-enhanced scoring; won CACHE Challenge #1; better enrichment than Vina | Open source |
| **DiffDock** | Diffusion-based blind docking; 38% top-1 success rate | Open source |
| **Uni-Dock** | GPU-accelerated; screened 38.2M compounds in 12 hours | Open source |
| **FEP-SPell-ABFE** | Free energy perturbation (gold standard accuracy) | Open source |

### 2.4 ADMET (M7) — Beyond Hand-Crafted Rules

| Tool | For | Improvement |
|------|-----|------------|
| **ADMETlab 3.0 API** (already coded!) | Small molecules | 119 endpoints vs. basic Lipinski; DMPNN architecture |
| **ToxinPred 3.0** | Peptide toxicity | ML model vs. GRAVY heuristic |
| **HemoPI 2.0** | Hemolysis prediction | ML model vs. charge+hydrophobicity rule |
| **CellPPD** | Cell penetration prediction | ML model vs. hand-crafted CPP score |
| **PeptideRanker** | Peptide bioactivity | Validated neural network |
| **NetMHCpan** | Immunogenicity | MHC binding prediction (missing entirely) |

### 2.5 Databases to Add

| Database | Purpose | Priority |
|----------|---------|----------|
| **BindingDB** | 3.2M quantitative binding data points | HIGH |
| **Open Targets** | Target-disease associations, 440M entities | HIGH |
| **PubChem** | Compound metadata, bioassay data | HIGH |
| **HCDT 2.0** | Validated drug-target interactions | MEDIUM |
| **Enamine REAL** | 1.4B make-on-demand compounds | MEDIUM |
| **CycPeptMPDB** | Cyclic peptide membrane permeability | MEDIUM |
| **STRING** | Protein-protein interaction networks | LOW |

---

## Part 3: Business Strategy

### 3.1 What Makes Platforms Commercially Viable

The top platforms (Insilico: $888M Servier deal; Iambic: $1.7B Takeda deal; Recursion: $1B+ valuation) differentiate on:

1. **Clinical-stage evidence** — Programs in Phase 1+ trials
2. **Speed** — IND nomination in <18 months (vs 3-4 years traditional)
3. **Proprietary data moats** — Unique training data beyond public databases
4. **Closed-loop automation** — Weekly design-make-test-learn cycles
5. **Validated hit rates** — >5x enrichment vs. random screening

### 3.2 Business Features Missing From Our Pipeline

| Feature | What It Provides | Business Impact |
|---------|-----------------|----------------|
| Patent landscape checking (SureChEMBL) | Freedom-to-operate analysis | Required for pharma partnerships |
| Synthesis accessibility scoring (SA Score) | Cost-of-goods estimation | Manufacturing feasibility |
| Retrosynthesis planning (ASKCOS/IBM RXN) | Synthetic route prediction | Practical candidate prioritization |
| PK/PD modeling | Dose prediction, therapeutic index | Clinical translatability |
| Experimental protocol generation | Specific assay recommendations | Faster wet-lab validation |
| Competitive landscape analysis | Candidates vs. drugs in clinical trials | Strategic positioning |
| Closed-loop feedback | Results from testing feed back to generation | Key commercial differentiator |

### 3.3 Revenue Models for AI Drug Discovery

| Model | Description | Our Fit |
|-------|-------------|---------|
| **Discovery partnerships** | Upfront + milestones + royalties with pharma | Best medium-term path |
| **Platform-as-a-Service** | Per-project fees for target analysis + candidate generation | Best short-term path |
| **Foundation model licensing** | License trained models to pharma ($50M+ deals) | Requires proprietary data |
| **Internal pipeline** | Own drug programs through clinical trials | Highest risk/reward |

### 3.4 Key Commercial Metrics to Track

| Metric | Target | Current Estimate |
|--------|--------|-----------------|
| Hit rate from AI screening | >5x enrichment | Unknown (no validation) |
| IND nomination time | <18 months | N/A |
| Molecules screened per IND | <500 | ~20-50 candidates generated |
| ADMET pass rate | >30% | Unknown |
| Synthesis success rate | >90% | Not assessed |

---

## Part 4: Prioritized Implementation Roadmap

### Phase 1: Critical Fixes (1-2 weeks)

- [ ] Add ADMET to composite scoring (weight 0.15)
- [ ] Remove `sequence_diversity` from composite (apply as post-ranking filter)
- [ ] Use `structure_confidence` as gate (>0.5) not weighted component
- [ ] Wire existing ADMETlab 3.0 client into `admet.py` for small molecules
- [ ] Add absolute binding energy threshold (discard > -5.0 kcal/mol)
- [ ] Re-dock modified candidates after M6

### Phase 2: Tool Upgrades (2-4 weeks)

- [ ] Integrate GNINA or DiffDock as docking alternative
- [ ] Integrate Boltz-1/Boltz-2 for structure prediction + affinity
- [ ] Replace peptide ADMET heuristics with ML models (ToxinPred, HemoPI, CellPPD)
- [ ] Add immunogenicity prediction (NetMHCpan)
- [ ] Add BindingDB and Open Targets databases
- [ ] Integrate RFpeptides for macrocyclic peptide design

### Phase 3: Scientific Rigor (1-2 months)

- [ ] Implement docking-based selectivity (dock against anti-target structures)
- [ ] Add MD simulation module (OpenMM) for top-N binding stability
- [ ] Add FEP calculations for lead compounds
- [ ] Implement real REINVENT 4 or SyntheMol for SM de novo design
- [ ] Build retrospective validation benchmarks
- [ ] Add MOA prediction capability

### Phase 4: Commercial Features (2-3 months)

- [ ] Patent landscape analysis (SureChEMBL integration)
- [ ] Synthetic accessibility scoring (RDKit SA Score)
- [ ] Retrosynthesis planning integration
- [ ] Experimental protocol generation
- [ ] Competitive landscape analysis
- [ ] PK/PD modeling integration
- [ ] Web API / dashboard for platform-as-a-service model

---

## Part 5: Estimated Candidate Reliability

### Current Pipeline

| Modality | Estimated In Vitro Hit Rate | Notes |
|----------|---------------------------|-------|
| Peptide | 5-15% (1-3 of 20 candidates) | PepMLM provides target-conditioned design; Vina docking + heuristic filtering |
| Small Molecule | ~0% without curation | Random fragment assembly; no target awareness |

### After Phase 2 Upgrades

| Modality | Projected Hit Rate | Improvement Driver |
|----------|-------------------|-------------------|
| Peptide | 20-35% | RFpeptides + GNINA + Boltz-2 affinity + ML ADMET |
| Small Molecule | 5-15% | Real generative model + GNINA + ADMETlab 3.0 |

### Industry Benchmark

| Platform | Reported Hit Rate |
|----------|-----------------|
| Insilico Medicine | 60-200 molecules to IND |
| Iambic Therapeutics | ~200 molecules to IND |
| Traditional HTS | 3,000-5,000 molecules to IND |

---

*Generated by strategic pipeline review — Brown Biotech Inc.*
