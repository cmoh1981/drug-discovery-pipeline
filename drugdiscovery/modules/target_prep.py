"""M1: Target Preparation.

Resolves a target identifier (gene name, UniProt ID, or PDB ID) to a
fully-populated TargetProfile, including sequence, structure, binding
pockets, anti-targets, and tissue classification.
"""

from __future__ import annotations

import dataclasses
import logging
import re
from pathlib import Path
from typing import Optional

from drugdiscovery.types import BindingPocket, PipelineConfig, TargetProfile
from drugdiscovery.utils.io import write_json
from drugdiscovery.utils.web import (
    fetch_alphafold_pdb,
    fetch_uniprot_fasta,
    fetch_uniprot_json,
    search_uniprot,
)

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Regex patterns for identifier classification
# ---------------------------------------------------------------------------

# UniProt accession: 6-10 alphanumeric characters starting with a letter.
# Matches both reviewed (P/Q/O prefix) and unreviewed forms.
_UNIPROT_RE = re.compile(r"^[A-Za-z][A-Za-z0-9]{5,9}$")

# PDB IDs are exactly 4 characters: digit + 3 alphanumeric.
_PDB_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")


# ---------------------------------------------------------------------------
# 1. resolve_target
# ---------------------------------------------------------------------------

def resolve_target(identifier: str) -> tuple[str, str]:
    """Resolve a target identifier to (gene_name, uniprot_id).

    Accepts:
      - UniProt accession  (e.g. "Q9Y2Z4")
      - PDB ID             (e.g. "1X8X")
      - Gene name          (e.g. "YARS2", "EGFR")

    Returns:
      (gene_name, uniprot_id) — both non-empty strings.

    Raises:
      ValueError if resolution fails.
    """
    ident = identifier.strip()

    # --- UniProt accession ---
    if _UNIPROT_RE.match(ident):
        logger.info("Identifier '%s' looks like a UniProt accession; fetching directly.", ident)
        data = fetch_uniprot_json(ident)
        gene_name = _extract_gene_name(data)
        logger.info("Resolved UniProt %s -> gene '%s'", ident, gene_name)
        return gene_name, ident

    # --- PDB ID ---
    if _PDB_RE.match(ident):
        logger.info("Identifier '%s' looks like a PDB ID; searching UniProt cross-ref.", ident)
        results = search_uniprot(f"xref_pdb:{ident}", limit=5)
        if not results:
            raise ValueError(f"No UniProt entry found for PDB ID '{ident}'")
        entry = results[0]
        uniprot_id = _entry_accession(entry)
        data = fetch_uniprot_json(uniprot_id)
        gene_name = _extract_gene_name(data)
        logger.info("Resolved PDB %s -> UniProt %s / gene '%s'", ident, uniprot_id, gene_name)
        return gene_name, uniprot_id

    # --- Gene name ---
    logger.info("Treating '%s' as a gene name; searching UniProt (human).", ident)
    results = search_uniprot(f"gene_exact:{ident}", limit=10)
    if not results:
        # Fallback to a broader search
        results = search_uniprot(ident, limit=10)
    if not results:
        raise ValueError(f"No UniProt entry found for gene name '{ident}'")

    # Prefer reviewed (Swiss-Prot) human entries
    entry = _pick_best_human_entry(results, ident)
    uniprot_id = _entry_accession(entry)
    data = fetch_uniprot_json(uniprot_id)
    gene_name = _extract_gene_name(data)
    logger.info("Resolved gene '%s' -> UniProt %s", ident, uniprot_id)
    return gene_name, uniprot_id


# ---------------------------------------------------------------------------
# 2. detect_binding_pockets
# ---------------------------------------------------------------------------

def detect_binding_pockets(uniprot_data: dict) -> list[BindingPocket]:
    """Parse UniProt JSON features to extract binding pockets.

    Groups nearby residues (within 5 positions) into the same pocket.

    Args:
        uniprot_data: Full UniProt JSON record.

    Returns:
        List of BindingPocket instances (may be empty).
    """
    target_feature_types = {"Active site", "Binding site", "Site"}
    features = uniprot_data.get("features", [])

    # Collect (residue_number, description, pocket_type) tuples
    raw: list[tuple[int, str, str]] = []
    for feat in features:
        ftype = feat.get("type", "")
        if ftype not in target_feature_types:
            continue
        loc = feat.get("location", {})
        start = loc.get("start", {}).get("value")
        end = loc.get("end", {}).get("value")
        if start is None:
            continue
        end = end if end is not None else start
        description = feat.get("description", "")
        pocket_type = _feature_type_to_pocket_type(ftype)
        for res in range(int(start), int(end) + 1):
            raw.append((res, description, pocket_type))

    if not raw:
        logger.info("No binding-site features found in UniProt data.")
        return []

    # Sort by residue number then group nearby residues into pockets
    raw.sort(key=lambda x: x[0])
    pockets: list[BindingPocket] = []
    _group_into_pockets(raw, pockets, gap_threshold=5)

    logger.info("Detected %d binding pocket(s) from UniProt features.", len(pockets))
    return pockets


# ---------------------------------------------------------------------------
# 3. find_anti_targets
# ---------------------------------------------------------------------------

def find_anti_targets(
    gene_name: str,
    uniprot_id: str,
    organism_id: str = "9606",
) -> list[dict]:
    """Find related proteins (same gene family) to use as anti-targets.

    Strips trailing digits from the gene name (e.g. "YARS2" -> "YARS") and
    searches for other family members in the same organism.

    Args:
        gene_name:   Primary target gene name.
        uniprot_id:  Primary target UniProt accession (excluded from results).
        organism_id: NCBI taxonomy ID (default 9606 = human).

    Returns:
        List of dicts with keys: gene_name, uniprot_id, similarity_note.
    """
    # Strip trailing digits to get the family prefix (YARS2 -> YARS, CDK1 -> CDK)
    family_prefix = re.sub(r"\d+$", "", gene_name).strip()
    if not family_prefix:
        family_prefix = gene_name

    logger.info(
        "Searching for anti-targets in family '%s' (excluding %s).",
        family_prefix, uniprot_id,
    )

    try:
        results = search_uniprot(
            f"gene:{family_prefix}*",
            limit=20,
            organism=organism_id,
        )
    except Exception as exc:
        logger.warning("Anti-target search failed: %s", exc)
        return []

    anti_targets: list[dict] = []
    for entry in results:
        acc = _entry_accession(entry)
        if acc == uniprot_id:
            continue
        # Extract gene name from the search result (may be abbreviated)
        entry_gene = _entry_gene_name(entry)
        similarity_note = (
            f"Same gene family as {gene_name} (prefix '{family_prefix}')"
        )
        anti_targets.append({
            "gene_name": entry_gene,
            "uniprot_id": acc,
            "similarity_note": similarity_note,
        })

    logger.info("Found %d anti-target candidate(s).", len(anti_targets))
    return anti_targets


# ---------------------------------------------------------------------------
# 4. determine_target_tissue
# ---------------------------------------------------------------------------

def determine_target_tissue(uniprot_data: dict) -> str:
    """Infer target tissue from subcellular location and GO terms.

    Args:
        uniprot_data: Full UniProt JSON record.

    Returns:
        A tissue string such as "systemic", "blood", "lung", "liver", etc.
    """
    location_text = _collect_location_text(uniprot_data).lower()

    # Priority-ordered keyword mapping
    mappings: list[tuple[str, str]] = [
        ("mitochondri",   "systemic"),
        ("nucleus",       "systemic"),
        ("cytoplasm",     "systemic"),
        ("cytosol",       "systemic"),
        ("endoplasmic",   "systemic"),
        ("golgi",         "systemic"),
        ("lysosom",       "systemic"),
        ("peroxisom",     "systemic"),
        ("extracellular", "blood"),
        ("secreted",      "blood"),
        ("plasma membrane", "systemic"),
        ("membrane",      "systemic"),
        ("lung",          "lung"),
        ("liver",         "liver"),
        ("kidney",        "kidney"),
        ("brain",         "brain"),
        ("heart",         "heart"),
        ("muscle",        "muscle"),
    ]

    for keyword, tissue in mappings:
        if keyword in location_text:
            logger.debug("Tissue '%s' inferred from keyword '%s'.", tissue, keyword)
            return tissue

    logger.info("Could not infer tissue from subcellular location; defaulting to 'systemic'.")
    return "systemic"


# ---------------------------------------------------------------------------
# 5. prepare_target (main entry point)
# ---------------------------------------------------------------------------

def prepare_target(cfg: PipelineConfig, output_dir: Path) -> TargetProfile:
    """M1 orchestrator: resolve target and build a full TargetProfile.

    Saves to output_dir:
      - target_profile.json
      - target.fasta
      - binding_pockets.json
      - anti_targets.json
      - structure.pdb  (if AlphaFold available)

    Args:
        cfg:        Pipeline configuration (must have cfg.target set).
        output_dir: Directory to write all M1 artefacts.

    Returns:
        Populated TargetProfile.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("M1 | Resolving target: '%s'", cfg.target)

    # Step 1: Resolve identifier -> (gene_name, uniprot_id)
    gene_name, uniprot_id = resolve_target(cfg.target)

    # Step 2: Fetch full UniProt JSON record
    logger.info("M1 | Fetching UniProt JSON for %s", uniprot_id)
    uniprot_data = fetch_uniprot_json(uniprot_id)

    # Step 3: Fetch FASTA sequence
    logger.info("M1 | Fetching FASTA sequence for %s", uniprot_id)
    fasta_text = fetch_uniprot_fasta(uniprot_id)
    sequence = _parse_sequence_from_fasta(fasta_text)

    # Step 4: Fetch structure (user-provided PDB or AlphaFold)
    structure_pdb_path: Optional[str] = None
    structure_source = ""
    if cfg.target_pdb_path and Path(cfg.target_pdb_path).is_file():
        pdb_path = output_dir / "structure.pdb"
        import shutil
        shutil.copy2(cfg.target_pdb_path, pdb_path)
        structure_pdb_path = str(pdb_path)
        structure_source = "user_pdb"
        logger.info("M1 | User-provided PDB copied from %s to %s", cfg.target_pdb_path, pdb_path)
    else:
        logger.info("M1 | Fetching AlphaFold structure for %s", uniprot_id)
        pdb_text = fetch_alphafold_pdb(uniprot_id)
        if pdb_text:
            pdb_path = output_dir / "structure.pdb"
            pdb_path.write_text(pdb_text, encoding="utf-8")
            structure_pdb_path = str(pdb_path)
            structure_source = "alphafold"
            logger.info("M1 | AlphaFold structure saved to %s", pdb_path)
        else:
            logger.warning("M1 | AlphaFold structure not available for %s; continuing.", uniprot_id)

    # Step 5: Detect binding pockets
    logger.info("M1 | Detecting binding pockets from UniProt features.")
    pockets = detect_binding_pockets(uniprot_data)

    # Step 6: Find anti-targets (optional)
    logger.info("M1 | Searching for anti-targets in same gene family.")
    try:
        anti_targets = find_anti_targets(gene_name, uniprot_id)
    except Exception as exc:
        logger.warning("M1 | Anti-target search failed (continuing): %s", exc)
        anti_targets = []

    # Step 7: Determine tissue
    if cfg.tissue:
        target_tissue = cfg.tissue
        logger.info("M1 | Using user-supplied tissue: '%s'", target_tissue)
    else:
        logger.info("M1 | Auto-detecting target tissue from UniProt data.")
        target_tissue = determine_target_tissue(uniprot_data)
        logger.info("M1 | Auto-detected tissue: '%s'", target_tissue)

    # Step 8: Extract additional metadata
    protein_name = _extract_protein_name(uniprot_data)
    function_description = _extract_function(uniprot_data)
    subcellular_location = _extract_subcellular_location_text(uniprot_data)
    go_terms = _extract_go_terms(uniprot_data)
    molecular_weight = _extract_molecular_weight(uniprot_data)

    # Step 9: Build TargetProfile
    profile = TargetProfile(
        gene_name=gene_name,
        uniprot_id=uniprot_id,
        organism="Homo sapiens",
        protein_name=protein_name,
        sequence=sequence,
        sequence_length=len(sequence),
        molecular_weight=molecular_weight,
        function_description=function_description,
        subcellular_location=subcellular_location,
        go_terms=go_terms,
        binding_pockets=pockets,
        structure_pdb_path=structure_pdb_path,
        structure_source=structure_source,
        anti_targets=anti_targets,
        target_tissue=target_tissue,
        metadata={
            "uniprot_query": cfg.target,
            "sequence_length": len(sequence),
            "pocket_count": len(pockets),
            "anti_target_count": len(anti_targets),
        },
    )

    # Step 10: Persist artefacts
    _save_artefacts(profile, fasta_text, output_dir)

    logger.info(
        "M1 | Target preparation complete: %s / %s (%d aa, %d pockets, %d anti-targets)",
        gene_name, uniprot_id, len(sequence), len(pockets), len(anti_targets),
    )
    return profile


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _extract_gene_name(uniprot_data: dict) -> str:
    """Pull the primary gene name from a UniProt JSON record."""
    try:
        return uniprot_data["genes"][0]["geneName"]["value"]
    except (KeyError, IndexError, TypeError):
        pass
    # Fallback: use the accession
    return uniprot_data.get("primaryAccession", "UNKNOWN")


def _extract_protein_name(uniprot_data: dict) -> str:
    """Pull the recommended full protein name."""
    try:
        return (
            uniprot_data["proteinDescription"]["recommendedName"]["fullName"]["value"]
        )
    except (KeyError, TypeError):
        return ""


def _extract_function(uniprot_data: dict) -> str:
    """Extract function text from comments."""
    for comment in uniprot_data.get("comments", []):
        if comment.get("commentType") == "FUNCTION":
            texts = comment.get("texts", [])
            if texts:
                return texts[0].get("value", "")
    return ""


def _extract_subcellular_location_text(uniprot_data: dict) -> str:
    """Build a human-readable subcellular location string."""
    parts: list[str] = []
    for comment in uniprot_data.get("comments", []):
        if comment.get("commentType") != "SUBCELLULAR LOCATION":
            continue
        for loc in comment.get("subcellularLocations", []):
            location_obj = loc.get("location", {})
            loc_value = location_obj.get("value", "")
            if loc_value:
                parts.append(loc_value)
    return "; ".join(parts) if parts else ""


def _extract_go_terms(uniprot_data: dict) -> list[str]:
    """Collect GO term descriptions from cross-references."""
    go_terms: list[str] = []
    for xref in uniprot_data.get("uniProtKBCrossReferences", []):
        if xref.get("database") != "GO":
            continue
        for prop in xref.get("properties", []):
            if prop.get("key") in ("GoTerm", "term"):
                go_terms.append(prop.get("value", ""))
    return go_terms


def _extract_molecular_weight(uniprot_data: dict) -> float:
    """Extract molecular weight in kDa (UniProt stores Da)."""
    try:
        mw_da = uniprot_data["sequence"]["molWeight"]
        return round(mw_da / 1000.0, 2)
    except (KeyError, TypeError):
        return 0.0


def _collect_location_text(uniprot_data: dict) -> str:
    """Aggregate all location-relevant text for tissue detection."""
    parts: list[str] = []

    # Subcellular location comments
    parts.append(_extract_subcellular_location_text(uniprot_data))

    # GO cellular_component terms (prefix C:)
    for term in _extract_go_terms(uniprot_data):
        if term.startswith("C:"):
            parts.append(term[2:])

    return " ".join(parts)


def _feature_type_to_pocket_type(feature_type: str) -> str:
    """Map UniProt feature type to our internal pocket_type vocabulary."""
    mapping = {
        "Active site": "active_site",
        "Binding site": "binding_site",
        "Site": "site",
    }
    return mapping.get(feature_type, "binding_site")


def _group_into_pockets(
    raw: list[tuple[int, str, str]],
    pockets: list[BindingPocket],
    gap_threshold: int = 5,
) -> None:
    """Greedy grouping: residues within gap_threshold of each other form one pocket."""
    if not raw:
        return

    # Build initial groups
    groups: list[list[tuple[int, str, str]]] = []
    current_group: list[tuple[int, str, str]] = [raw[0]]

    for item in raw[1:]:
        resnum = item[0]
        last_resnum = current_group[-1][0]
        if resnum - last_resnum <= gap_threshold:
            current_group.append(item)
        else:
            groups.append(current_group)
            current_group = [item]
    groups.append(current_group)

    for idx, group in enumerate(groups, start=1):
        residue_numbers = sorted({r for r, _, _ in group})
        # Use the most common description in the group
        desc_counter: dict[str, int] = {}
        type_counter: dict[str, int] = {}
        for _, desc, ptype in group:
            desc_counter[desc] = desc_counter.get(desc, 0) + 1
            type_counter[ptype] = type_counter.get(ptype, 0) + 1

        description = max(desc_counter, key=lambda k: desc_counter[k])
        pocket_type = max(type_counter, key=lambda k: type_counter[k])

        pockets.append(
            BindingPocket(
                pocket_id=f"pocket_{idx:02d}",
                residue_numbers=residue_numbers,
                description=description,
                pocket_type=pocket_type,
            )
        )


def _parse_sequence_from_fasta(fasta_text: str) -> str:
    """Extract bare amino acid sequence from raw FASTA text."""
    lines = fasta_text.strip().splitlines()
    seq_parts: list[str] = []
    for line in lines:
        if line.startswith(">"):
            continue
        seq_parts.append(line.strip())
    return "".join(seq_parts)


def _entry_accession(entry: dict) -> str:
    """Extract accession from a UniProt search result entry."""
    # Search API returns primaryAccession
    acc = entry.get("primaryAccession", "")
    if acc:
        return acc
    # Fallback for alternative response shapes
    return entry.get("accession", "")


def _entry_gene_name(entry: dict) -> str:
    """Extract gene name from a UniProt search result entry (best-effort)."""
    try:
        return entry["genes"][0]["geneName"]["value"]
    except (KeyError, IndexError, TypeError):
        pass
    # gene_names field from search API (space-separated)
    raw = entry.get("gene_names") or entry.get("geneNames", "")
    if isinstance(raw, str) and raw:
        return raw.split()[0]
    return entry.get("primaryAccession", "UNKNOWN")


def _pick_best_human_entry(results: list[dict], gene_name: str) -> dict:
    """From search results, prefer the entry whose gene name matches exactly."""
    gene_upper = gene_name.upper()
    for entry in results:
        candidate_gene = _entry_gene_name(entry).upper()
        if candidate_gene == gene_upper:
            return entry
    # If no exact match, return the first result
    return results[0]


def _save_artefacts(
    profile: TargetProfile,
    fasta_text: str,
    output_dir: Path,
) -> None:
    """Persist all M1 output files."""
    # target_profile.json — convert dataclass to dict, handle nested dataclasses
    profile_dict = dataclasses.asdict(profile)
    write_json(output_dir / "target_profile.json", profile_dict)

    # FASTA
    fasta_path = output_dir / "target.fasta"
    fasta_path.write_text(fasta_text, encoding="utf-8")
    logger.info("Saved FASTA to %s", fasta_path)

    # binding_pockets.json
    pockets_data = [dataclasses.asdict(p) for p in profile.binding_pockets]
    write_json(output_dir / "binding_pockets.json", pockets_data)

    # anti_targets.json
    write_json(output_dir / "anti_targets.json", profile.anti_targets)

    logger.info("M1 | All artefacts saved to %s", output_dir)
