"""Core data types for the drug discovery pipeline."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional


class Modality(str, Enum):
    SMALL_MOLECULE = "small_molecule"
    PEPTIDE = "peptide"


class ModeOfAction(str, Enum):
    AGONIST = "agonist"
    ANTAGONIST = "antagonist"


@dataclass
class BindingPocket:
    """A detected binding pocket on the target protein."""
    pocket_id: str
    residue_numbers: list[int]
    residue_names: list[str] = field(default_factory=list)
    description: str = ""
    pocket_type: str = ""  # e.g. "active_site", "binding_site", "allosteric"
    centroid: Optional[tuple[float, float, float]] = None


@dataclass
class TargetProfile:
    """Consolidated target protein information."""
    gene_name: str
    uniprot_id: str
    organism: str = "Homo sapiens"
    protein_name: str = ""
    sequence: str = ""
    sequence_length: int = 0
    molecular_weight: float = 0.0
    function_description: str = ""
    subcellular_location: str = ""
    go_terms: list[str] = field(default_factory=list)
    binding_pockets: list[BindingPocket] = field(default_factory=list)
    structure_pdb_path: Optional[str] = None
    structure_source: str = ""  # "alphafold", "pdb", "esmfold"
    anti_targets: list[dict] = field(default_factory=list)
    target_tissue: str = ""
    metadata: dict = field(default_factory=dict)


@dataclass
class Candidate:
    """Universal candidate schema for both small molecules and peptides."""
    candidate_id: str
    candidate_type: str = ""      # "library_hit", "de_novo", "modified"
    modality: str = ""            # "small_molecule", "peptide"
    source: str = ""              # "chembl", "pepmlm", "reinvent4", etc.
    smiles: str = ""              # SMILES (SM only)
    iupac_name: str = ""          # IUPAC name (SM only, resolved via PubChem)
    sequence: str = ""            # AA sequence (peptide only)
    molecular_weight: float = 0.0
    net_charge: float = 0.0
    binding_score: float = 0.0    # Best binding score (kcal/mol or normalized)
    moa_predicted: str = "unknown"  # "agonist", "antagonist", "unknown"
    selectivity_score: float = 0.0  # 0-1 vs anti-targets
    drug_likeness: float = 0.0     # 0-1
    admet_score: float = 0.0       # 0-1 aggregate
    admet_flags: int = 0           # Count of ADMET issues
    delivery_system: str = ""      # Recommended delivery
    composite_score: float = 0.0   # 0-1 weighted final score
    rank: int = 0                  # Final rank (1 = best)
    # Additional properties
    gravy: float = 0.0
    isoelectric_point: float = 0.0
    perplexity: float = 0.0
    structure_confidence: float = 0.0
    sequence_diversity: float = 0.0
    parent_id: str = ""            # For modified candidates
    modification: str = ""         # Modification type applied
    modification_detail: str = ""
    # Synthetic accessibility & pose validation
    sa_score: float = 0.0                 # RDKit SA score (1-10, lower = easier to synthesize)
    posecheck_valid: bool = True          # PoseCheck validation result
    # Perturbation biology (M4.6)
    perturbation_score: float = 0.0       # 0-1 aggregate perturbation score
    cmap_connectivity: float = 0.0        # L1000/CMAP connectivity score (-1 to 1, normalized to 0-1)
    cmap_compound_match: str = ""         # Nearest CMAP compound matched
    network_effect_score: float = 0.0     # STRING PPI network propagation score (0-1)
    disease_signature_reversal: float = 0.0  # How well drug reverses disease signature (0-1)
    metadata: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        """Convert to flat dictionary for CSV export."""
        d = {
            "candidate_id": self.candidate_id,
            "candidate_type": self.candidate_type,
            "modality": self.modality,
            "source": self.source,
            "smiles": self.smiles,
            "iupac_name": self.iupac_name,
            "sequence": self.sequence,
            "molecular_weight": round(self.molecular_weight, 2),
            "net_charge": round(self.net_charge, 2),
            "binding_score": round(self.binding_score, 4),
            "moa_predicted": self.moa_predicted,
            "selectivity_score": round(self.selectivity_score, 4),
            "drug_likeness": round(self.drug_likeness, 4),
            "admet_score": round(self.admet_score, 4),
            "admet_flags": self.admet_flags,
            "delivery_system": self.delivery_system,
            "perturbation_score": round(self.perturbation_score, 4),
            "cmap_connectivity": round(self.cmap_connectivity, 4),
            "cmap_compound_match": self.cmap_compound_match,
            "network_effect_score": round(self.network_effect_score, 4),
            "disease_signature_reversal": round(self.disease_signature_reversal, 4),
            "sa_score": round(self.sa_score, 2),
            "composite_score": round(self.composite_score, 4),
            "rank": self.rank,
        }
        return d


@dataclass
class ADMETProfile:
    """ADMET prediction results for a candidate."""
    candidate_id: str
    # Absorption
    solubility: float = 0.0          # 0-1
    permeability: float = 0.0       # 0-1
    cpp_score: float = 0.0          # Cell-penetrating peptide score (peptide)
    oral_bioavailability: float = 0.0  # (SM)
    # Distribution
    bbb_permeability: float = 0.0
    plasma_protein_binding: float = 0.0
    # Metabolism
    protease_stability: float = 0.0  # Peptide
    cyp_inhibition: float = 0.0      # SM
    half_life_estimate: str = ""
    # Excretion
    renal_clearance: str = ""
    # Toxicity
    hemolysis_risk: float = 0.0
    aggregation_propensity: float = 0.0
    hepatotoxicity_risk: float = 0.0
    herg_liability: float = 0.0      # SM cardiac risk
    # Aggregate
    flag_count: int = 0
    aggregate_score: float = 0.0     # 0-1
    flags: list[str] = field(default_factory=list)


@dataclass
class DeliveryRecommendation:
    """Delivery system recommendation for a candidate."""
    candidate_id: str
    primary_system: str = ""
    secondary_system: str = ""
    rationale: str = ""
    tissue: str = ""
    route: str = ""           # "oral", "IV", "subcutaneous", etc.
    formulation_notes: str = ""


@dataclass
class PipelineConfig:
    """Pipeline run configuration."""
    modality: Modality = Modality.PEPTIDE
    mode: ModeOfAction = ModeOfAction.ANTAGONIST
    target: str = ""
    tissue: str = ""
    top_n: int = 20
    device: str = "auto"
    use_runpod: bool = False
    output_dir: str = "results/"
    config_file: str = ""
    target_pdb_path: str = ""  # User-provided PDB file path
    # Module-specific settings loaded from YAML
    scoring_weights: dict = field(default_factory=lambda: {
        "binding_energy": 0.30,
        "selectivity": 0.20,
        "drug_likeness": 0.10,
        "admet_aggregate": 0.15,
        "moa_consistency": 0.10,
        "perturbation": 0.15,
    })
    pepmlm_settings: dict = field(default_factory=lambda: {
        "model_name": "TianlaiChen/PepMLM-650M",
        "lengths": [8, 10, 12, 15, 20, 25, 30],
        "num_per_length": 50,
        "top_k": 3,
    })
    diffpepbuilder_settings: dict = field(default_factory=lambda: {
        "num_samples": 100,
        "peptide_length_min": 8,
        "peptide_length_max": 25,
    })
    vina_settings: dict = field(default_factory=lambda: {
        "exhaustiveness": 32,
        "n_poses": 10,
        "box_size": 30.0,
    })
    library_settings: dict = field(default_factory=lambda: {
        "max_hits_per_db": 500,
        "similarity_threshold": 0.7,
    })
