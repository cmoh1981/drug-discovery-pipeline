# Updates — 2026-02-27

## Completed Today

### 1. Vina Docking Integration (M4.5)

Wired AutoDock Vina molecular docking into the pipeline as a new step between Structure Prediction (M4) and Binding & Scoring (M5). Previously, all `binding_score` values were 0.0 — the 30% binding_energy weight in the composite score contributed nothing.

**Files changed:**

| File | Change |
|------|--------|
| `drugdiscovery/modules/docking.py` | **NEW** — `dock_candidates()` orchestrates Vina docking for all candidates |
| `drugdiscovery/pipeline.py` | Inserted M4.5 step between M4 and M5 |
| `drugdiscovery/types.py` | Added `box_size: 30.0` to `vina_settings` default |
| `requirements.txt` | Uncommented `vina>=1.2.0` and `meeko>=0.5.0` |

**How it works:**
- Checks for a receptor PDB (`target_profile.structure_pdb_path`)
- Determines docking center: pocket centroid > pocket residue CA average > whole-protein CA average
- Small molecules: calls `tools.vina.dock_smiles()`
- Peptides: calls `tools.vina.dock_peptide_pdb()` using conformer/structure PDB from metadata
- Writes `docking/docking_results.csv` to run directory
- Graceful skip if no PDB, no vina/meeko, or no docking center available

### 2. ESRRG Results Document

Created `ESRRG.md` at project root summarizing the first pipeline run (ESRRG small-molecule agonist, run ID `20260227_134323`). Includes target profile, top 20 rankings, ADMET summary, key observations, and next steps.

### 3. Skill File Updated

Updated `~/.claude/skills/drugdiscovery/SKILL.md` to document the 10-step pipeline (M1–M9 + M4.5), docking graceful degradation, and added meeko/vina to status checks. Future Claude sessions will know about M4.5.

---

## Current Pipeline State

```
M1  Target Preparation        ✅ Working
M2  Library Screening          ✅ Working (ChEMBL only — see issues below)
M3  De Novo Design             ✅ Working (RDKit enumeration — REINVENT4 not installed)
M4  Structure Prediction       ✅ Working (RDKit ETKDG conformers)
M4.5 Molecular Docking         ✅ NEW — implemented, untested with real Vina run
M5  Binding & Scoring          ✅ Working
M6  Modifications              ✅ Working (prodrug + salt forms for SM)
M7  ADMET Prediction           ✅ Working (rule-based)
M8  Delivery System            ✅ Working
M9  Reporting                  ✅ Working
```

**Modules:** 13 files in `drugdiscovery/modules/`
**Tools:** 6 tool wrappers in `drugdiscovery/tools/` (vina, pepmlm, colabfold, diffpepbuilder, deepdta, admetlab3)
**Databases:** 5 database connectors (chembl, drugbank, coconut, zinc22, peptide_dbs)

---

## Known Issues

### Database Access Failures (from ESRRG run log)

| Database | Error | Impact |
|----------|-------|--------|
| DrugBank | 403 Forbidden | No hits — needs API key or local CSV download |
| COCONUT | 404 Not Found | No hits — API endpoint may have changed |
| ZINC22 | Read timeout + JSON parse error | No hits — server unresponsive |

Only **ChEMBL** returned results (327 hits, 253 after dedup). Library coverage is limited to one database.

### No Binding Pocket Detected

UniProt returned zero binding site annotations for ESRRG. PDB 1S9P has a known ligand-binding domain (LBD) cavity, but we need pocket detection tools to find it automatically.

### De Novo Generation Is Fragment-Level

RDKit fragment enumeration produces very small molecules (152–286 Da). These are useful starting points but lack the complexity for potent binding. REINVENT4 would generate more drug-like candidates.

---

## Next Work

### Priority 1: Validate Docking (M4.5)

- [ ] Install `vina` and `meeko` packages: `pip install vina meeko`
- [ ] Re-run ESRRG pipeline: `/drugdiscovery run small_molecule agonist ESRRG --config custom_with_pdb.yaml`
- [ ] Confirm `binding_score` values are non-zero in output
- [ ] Verify `docking/docking_results.csv` exists with real scores
- [ ] Check that composite scores shift with real binding data (30% weight)

### Priority 2: Add Pocket Detection

- [ ] Integrate fpocket or P2Rank to auto-detect binding pockets from PDB
- [ ] Create `drugdiscovery/tools/pocket_detection.py`
- [ ] Wire into M1 (Target Preparation) to populate `binding_pockets` with centroid coordinates
- [ ] Test on PDB 1S9P — should find the ESRRG LBD cavity
- [ ] Eliminates fallback to whole-protein centroid for docking

### Priority 3: Fix Database Access

- [ ] **DrugBank**: Download open-data CSV from https://go.drugbank.com/releases/latest, set `DRUGBANK_DATA_DIR` env var
- [ ] **COCONUT**: Investigate new API endpoint (old `/api/search/simple` returns 404)
- [ ] **ZINC22**: Add longer timeout or switch to ZINC20 REST API as fallback
- [ ] Goal: at least 3/4 databases returning hits

### Priority 4: Improve De Novo Generation

- [ ] Install REINVENT4 for generative chemistry (replaces RDKit fragment enumeration)
- [ ] Or integrate MolGPT / other generative model for SM candidates
- [ ] Target: molecules in 300–500 Da range with drug-like scaffolds

### Priority 5: Enhance ADMET

- [ ] Wire in ADMETlab3 API (`tools/admetlab3.py` exists but may not be called)
- [ ] Replace or supplement rule-based ADMET with ML predictions
- [ ] Add CYP450 interaction predictions for SM candidates

### Priority 6: Peptide Pipeline Validation

- [ ] Run a peptide pipeline end-to-end (e.g., `/drugdiscovery run peptide antagonist YARS2`)
- [ ] Verify PepMLM de novo generation works with GPU
- [ ] Test peptide docking via `dock_peptide_pdb()`
- [ ] Check ColabFold/ESMFold structure prediction

---

## Run History

| Date | Target | Modality | Mode | Run ID | Notes |
|------|--------|----------|------|--------|-------|
| 2026-02-27 | ESRRG | small_molecule | agonist | `20260227_132851` | Early run (partial) |
| 2026-02-27 | ESRRG | small_molecule | agonist | `20260227_134323` | Full run, binding_score=0.0 (pre-docking) |

---

## File Tree (key files)

```
drug-discovery-pipeline/
  ESRRG.md                          # ESRRG run results summary
  updates0227.md                    # This file
  requirements.txt                  # Updated: vina + meeko uncommented
  drugdiscovery/
    pipeline.py                     # Updated: M4.5 docking step added
    types.py                        # Updated: box_size in vina_settings
    modules/
      docking.py                    # NEW: M4.5 Vina docking orchestration
      target_prep.py                # M1
      library_screening.py          # M2
      de_novo_sm.py                 # M3 (small molecule)
      de_novo_peptide.py            # M3 (peptide)
      structure_pred.py             # M4
      scoring.py                    # M5
      modifications.py              # M6
      admet.py                      # M7
      delivery.py                   # M8
      reporting.py                  # M9
    tools/
      vina.py                       # dock_smiles(), dock_peptide_pdb()
      pepmlm.py                     # PepMLM de novo peptide generation
      colabfold.py                  # ColabFold structure prediction
      diffpepbuilder.py             # DiffPepBuilder peptide design
      deepdta.py                    # DeepDTA binding affinity
      admetlab3.py                  # ADMETlab3 API client
    databases/
      chembl.py                     # ChEMBL (working)
      drugbank.py                   # DrugBank (403 error)
      coconut.py                    # COCONUT (404 error)
      zinc22.py                     # ZINC22 (timeout)
      peptide_dbs.py                # Peptide databases
  results/
    20260227_134323_ESRRG_small_molecule_agonist/   # Latest run
```
