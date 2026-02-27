# ESRRG Small-Molecule Agonist Discovery Results

**Run Date:** 2026-02-27
**Pipeline Version:** 0.1.0
**Run ID:** `20260227_134323_ESRRG_small_molecule_agonist`
**Elapsed:** 133.4 seconds

---

## Target Profile

| Property | Value |
|----------|-------|
| Gene | ESRRG |
| Protein | Estrogen-related receptor gamma |
| UniProt | [P62508](https://www.uniprot.org/uniprot/P62508) |
| Organism | Homo sapiens |
| Length | 458 aa (51.3 kDa) |
| Localization | Nucleus |
| Structure | PDB 1S9P (user-provided) |
| Tissue | Systemic |

**Function:** Orphan nuclear receptor that acts as a constitutive transcription activator. Binds estrogen response elements (EREs) and activates transcription without a bound ligand. Induces PERM1 expression in skeletal muscle and plays a role in cold-induced thermogenesis.

**GO Annotations (selected):**
- DNA-binding transcription activator (RNA Pol II)
- Estrogen response element binding
- Nuclear steroid receptor activity
- Positive regulation of cold-induced thermogenesis

**Anti-targets identified:** 4 isoforms (C9JU32, C9J0E3, C9J5W9, C9JNX5)
**Binding pockets detected:** 0 (UniProt features returned none)

---

## Pipeline Configuration

| Parameter | Value |
|-----------|-------|
| Modality | Small molecule |
| Mode of action | Agonist |
| Device | CPU |
| Top N | 20 |

### Scoring Weights

| Component | Weight |
|-----------|--------|
| Binding energy | 0.30 |
| Structure confidence | 0.20 |
| Selectivity | 0.20 |
| Drug-likeness | 0.15 |
| Sequence diversity | 0.15 |

---

## Pipeline Summary

| Stage | Description | Result |
|-------|-------------|--------|
| M1 | Target Preparation | ESRRG resolved via UniProt; PDB 1S9P loaded |
| M2 | Library Screening | 253 unique hits (327 ChEMBL, 0 DrugBank, 0 COCONUT, 0 ZINC22) |
| M3 | De Novo Design | 150 molecules (RDKit fragment enumeration; REINVENT4 unavailable) |
| M4 | Structure Prediction | 403/403 SM conformers generated (RDKit ETKDG) |
| M5 | Binding & Scoring | 403 candidates scored; binding_score range [0.0, 0.0] |
| M6 | Modifications | 40 variants (prodrug + salt forms of top 20) |
| M7 | ADMET Prediction | 60 candidates profiled |
| M8 | Delivery | All candidates: Oral delivery (systemic target) |
| M9 | Reporting | Final report generated |

**Total candidates evaluated:** 403 (initial) -> 60 (top 20 + 40 modifications)

> **Note:** Binding scores are all 0.0 because Vina docking was not wired into the pipeline at run time. The newly added M4.5 docking module will populate real binding affinities on the next run.

---

## Top 10 Candidates

### Rank 1: CHEMBL1671971 (Library Hit)

| Property | Value |
|----------|-------|
| Source | ChEMBL |
| SMILES | `COc1cc(/C=C2\SC(=O)NC2=O)ccc1Oc1ccc(C#N)cc1C(F)(F)F` |
| MW | 420.0 Da |
| Composite Score | 0.5367 |
| Drug-likeness | 1.00 (Lipinski Ro5: 0 violations) |
| ADMET Score | 0.62 (2 flags: hERG risk, hepatotoxicity risk) |
| Selectivity | 0.50 |
| Delivery | Oral |

This thiazolidinedione-like compound with a trifluoromethyl-cyanophenyl ether moiety ranked first. Known ChEMBL entry with activity data against ESRRG. Two ADMET flags (hERG liability and hepatotoxicity risk) warrant medicinal chemistry optimization.

**Variants generated:**
- CHEMBL1671971_PRODRUG1 (Rank 2) -- Ester prodrug strategy
- CHEMBL1671971_SALT1 (Rank 3) -- HCl salt form

### Rank 4: RDKIT_00143 (De Novo)

| Property | Value |
|----------|-------|
| Source | RDKit enumeration |
| SMILES | `NS(=O)(=O)N1CCNCC1` |
| MW | 165.1 Da |
| Composite Score | 0.5296 |
| Drug-likeness | 1.00 |
| ADMET Score | 0.78 (0 flags) |
| Selectivity | 0.50 |
| Delivery | Oral |

Piperazine sulfonamide fragment. Clean ADMET profile, low MW (fragment-like). Good starting point for fragment-based optimization.

### Rank 7: RDKIT_00125 (De Novo)

| Property | Value |
|----------|-------|
| Source | RDKit enumeration |
| SMILES | `FC(F)(F)COC1OCCC1SN1CCNCC1` |
| MW | 286.1 Da |
| Composite Score | 0.5231 |
| Drug-likeness | 1.00 |
| ADMET Score | 0.74 (0 flags) |
| BBB Permeability | 0.70 |
| Delivery | Oral |

Trifluoromethyl-containing heterocyclic with piperazine. Good BBB permeability, clean safety profile.

### Rank 10: RDKIT_00118 (De Novo)

| Property | Value |
|----------|-------|
| Source | RDKit enumeration |
| SMILES | `C1CCC(SC2CNCCN2)CC1` |
| MW | 200.1 Da |
| Composite Score | 0.5227 |
| Drug-likeness | 1.00 |
| ADMET Score | 0.74 (0 flags) |
| Delivery | Oral |

Cyclohexyl-thio-piperazine scaffold. Fragment-like MW, clean profile, BBB-permeable.

---

## Full Top 20 Rankings

| Rank | Candidate ID | Type | Source | MW (Da) | Score | ADMET | Flags |
|------|-------------|------|--------|---------|-------|-------|-------|
| 1 | CHEMBL1671971 | library_hit | chembl | 420.0 | 0.5367 | 0.62 | 2 |
| 2 | CHEMBL1671971_PRODRUG1 | modified | chembl | 420.0 | 0.5367 | 0.62 | 2 |
| 3 | CHEMBL1671971_SALT1 | modified | chembl | 420.0 | 0.5367 | 0.62 | 2 |
| 4 | RDKIT_00143 | de_novo | rdkit_enum | 165.1 | 0.5296 | 0.78 | 0 |
| 5 | RDKIT_00143_PRODRUG1 | modified | rdkit_enum | 165.1 | 0.5296 | 0.78 | 0 |
| 6 | RDKIT_00143_SALT1 | modified | rdkit_enum | 165.1 | 0.5296 | 0.78 | 0 |
| 7 | RDKIT_00125 | de_novo | rdkit_enum | 286.1 | 0.5231 | 0.74 | 0 |
| 8 | RDKIT_00125_PRODRUG1 | modified | rdkit_enum | 286.1 | 0.5231 | 0.74 | 0 |
| 9 | RDKIT_00125_SALT1 | modified | rdkit_enum | 286.1 | 0.5231 | 0.74 | 0 |
| 10 | RDKIT_00118 | de_novo | rdkit_enum | 200.1 | 0.5227 | 0.74 | 0 |
| 11 | RDKIT_00118_PRODRUG1 | modified | rdkit_enum | 200.1 | 0.5227 | 0.74 | 0 |
| 12 | RDKIT_00118_SALT1 | modified | rdkit_enum | 200.1 | 0.5227 | 0.74 | 0 |
| 13 | RDKIT_00074 | de_novo | rdkit_enum | 213.2 | 0.5224 | 0.74 | 0 |
| 14 | RDKIT_00074_PRODRUG1 | modified | rdkit_enum | 213.2 | 0.5224 | 0.74 | 0 |
| 15 | RDKIT_00074_SALT1 | modified | rdkit_enum | 213.2 | 0.5224 | 0.74 | 0 |
| 16 | RDKIT_00041 | de_novo | rdkit_enum | 228.2 | 0.5218 | 0.68 | 0 |
| 17 | RDKIT_00041_PRODRUG1 | modified | rdkit_enum | 228.2 | 0.5218 | 0.68 | 0 |
| 18 | RDKIT_00041_SALT1 | modified | rdkit_enum | 228.2 | 0.5218 | 0.68 | 0 |
| 19 | RDKIT_00033 | de_novo | rdkit_enum | 221.2 | 0.5216 | 0.74 | 0 |
| 20 | RDKIT_00033_PRODRUG1 | modified | rdkit_enum | 221.2 | 0.5216 | 0.74 | 0 |

---

## ADMET Summary

| Metric | Value |
|--------|-------|
| Mean ADMET score | 0.733 |
| Candidates with 0 flags | 57/60 (95%) |
| Candidates with flags | 3/60 (CHEMBL1671971 + variants) |
| Flagged risks | hERG liability, hepatotoxicity |
| All candidates | Oral delivery route |

---

## Key Observations

1. **Binding scores are unvalidated.** All binding_score values are 0.0 because docking was not yet integrated when this run executed. The 30% binding_energy weight in the composite score contributed nothing, meaning rankings were driven by drug-likeness (15%), selectivity (20%), structure confidence (20%), and sequence diversity (15%).

2. **CHEMBL1671971 is the only library hit in the top 20.** This thiazolidinedione derivative with a cyanophenyl ether is a known ESRRG-active compound from ChEMBL. All other top candidates are de novo RDKit enumeration products or their modifications.

3. **De novo candidates are fragment-like.** MW range 152-286 Da, well below 500 Da. These are suitable starting points for fragment growing/merging strategies but may lack sufficient complexity for potent binding.

4. **ADMET profiles are generally clean.** 95% of candidates have zero ADMET flags. Only CHEMBL1671971 (and its variants) show hERG and hepatotoxicity risks.

5. **No binding pocket was detected.** UniProt returned no binding site annotations for ESRRG. Future runs should use pocket detection tools (fpocket, P2Rank) on PDB 1S9P to identify the ligand-binding domain (LBD) cavity for targeted docking.

6. **Database coverage was limited.** Only ChEMBL returned hits (327). DrugBank (403 error), COCONUT (404 error), and ZINC22 (timeout) all failed. Improved database access would expand the chemical search space.

---

## Next Steps

1. **Re-run with docking enabled.** The M4.5 docking module is now integrated. Re-run with the same config to get real Vina binding affinities against PDB 1S9P.

2. **Add pocket detection.** Run fpocket or P2Rank on 1S9P to identify the LBD pocket and its centroid coordinates for focused docking.

3. **Expand de novo generation.** Install REINVENT4 for more sophisticated generative chemistry instead of RDKit fragment enumeration.

4. **Resolve database access.** Fix DrugBank authentication and COCONUT/ZINC22 endpoints to broaden library screening.

5. **Prioritize CHEMBL1671971 for optimization.** Address hERG and hepatotoxicity flags through medicinal chemistry (reduce lipophilicity, remove/modify the CF3 group).

---

## Output Files

```
results/20260227_134323_ESRRG_small_molecule_agonist/
  config_snapshot.yaml          # Run configuration
  pipeline.log                  # Full execution log
  target/
    target_profile.json         # ESRRG target data
    target.fasta                # Protein sequence
    structure.pdb               # PDB 1S9P
    binding_pockets.json        # (empty - none detected)
    anti_targets.json           # 4 isoforms
  library_hits/
    chembl_hits.csv             # 327 ChEMBL hits
    merged_hits.csv             # 253 deduplicated
  de_novo/
    candidates.csv              # 150 RDKit-generated molecules
    conformers/                 # 403 PDB conformer files
  scoring/
    scored_candidates.csv       # All 60 final candidates
    top_20_candidates.csv       # Top 20 ranked
  modifications/
    modified_candidates.csv     # 40 prodrug/salt variants
  admet/
    admet_profiles.csv          # ADMET predictions
  delivery/
    delivery_recommendations.csv
  report/
    final_report.md             # Detailed pipeline report
    executive_summary.csv       # Full ranked list
    figures/                    # Visualization plots
```
