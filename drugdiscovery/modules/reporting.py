"""M9: Final report generation.

Produces HTML, Markdown, and CSV reports with visualizations.
Refactored from YARS2 pipeline: analyze_results.py
"""

from __future__ import annotations

import html
import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Optional

from drugdiscovery.types import Candidate, PipelineConfig, TargetProfile
from drugdiscovery.tools.pubchem import resolve_iupac_names
from drugdiscovery.utils.io import write_csv

logger = logging.getLogger(__name__)


def generate_report(
    cfg: PipelineConfig,
    candidates: list[Candidate],
    target_profile: TargetProfile,
    output_dir: Path,
) -> None:
    """Generate final pipeline report in multiple formats."""
    output_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = output_dir / "figures"
    fig_dir.mkdir(exist_ok=True)

    # Resolve IUPAC names for small-molecule candidates
    resolve_iupac_names(candidates)

    # Executive summary CSV
    summary_rows = [c.to_dict() for c in candidates]
    write_csv(output_dir / "executive_summary.csv", summary_rows)

    # Generate visualizations
    _generate_figures(candidates, fig_dir)

    # Markdown report
    md_text = _build_markdown_report(cfg, candidates, target_profile)
    (output_dir / "final_report.md").write_text(md_text, encoding="utf-8")

    # HTML report
    html_text = _build_html_report(cfg, candidates, target_profile)
    (output_dir / "final_report.html").write_text(html_text, encoding="utf-8")

    logger.info("Reports saved to %s", output_dir)


# ---------------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------------

def _build_markdown_report(
    cfg: PipelineConfig,
    candidates: list[Candidate],
    target: TargetProfile,
) -> str:
    """Build full Markdown report."""
    lines = [
        f"# Drug Discovery Pipeline Report",
        f"",
        f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"",
        f"## Target Profile",
        f"",
        f"| Property | Value |",
        f"|----------|-------|",
        f"| Gene | {target.gene_name} |",
        f"| UniProt | {target.uniprot_id} |",
        f"| Organism | {target.organism} |",
        f"| Protein | {target.protein_name} |",
        f"| Length | {target.sequence_length} aa |",
        f"| Tissue | {target.target_tissue} |",
        f"| Binding Pockets | {len(target.binding_pockets)} |",
        f"| Anti-targets | {len(target.anti_targets)} |",
        f"",
        f"## Pipeline Configuration",
        f"",
        f"| Parameter | Value |",
        f"|-----------|-------|",
        f"| Modality | {cfg.modality.value} |",
        f"| Mode of Action | {cfg.mode.value} |",
        f"| Device | {cfg.device} |",
        f"| Top N | {cfg.top_n} |",
        f"",
        f"## Scoring Methodology",
        f"",
        f"| Component | Weight |",
        f"|-----------|--------|",
    ]
    for comp, wt in cfg.scoring_weights.items():
        lines.append(f"| {comp} | {wt:.2f} |")

    lines.extend([
        f"",
        f"## Top Candidates",
        f"",
    ])

    # Top candidates table
    top = candidates[:cfg.top_n]
    if top:
        header = "| Rank | ID | Name | Type | Source | Score | Drug-likeness | ADMET | Delivery |"
        sep = "|------|-----|------|------|--------|-------|---------------|-------|----------|"
        lines.append(header)
        lines.append(sep)
        for c in top:
            display_name = c.iupac_name[:30] if c.iupac_name else (c.sequence[:20] if c.sequence else "")
            lines.append(
                f"| {c.rank} | {c.candidate_id} | {display_name} | {c.candidate_type} | {c.source} "
                f"| {c.composite_score:.4f} | {c.drug_likeness:.2f} "
                f"| {c.admet_score:.2f} | {c.delivery_system} |"
            )

    # Detailed profiles
    lines.extend([
        f"",
        f"## Detailed Candidate Profiles",
        f"",
    ])
    for c in top[:10]:
        lines.extend([
            f"### {c.candidate_id} (Rank {c.rank})",
            f"",
            f"- **Type:** {c.candidate_type} ({c.source})",
            f"- **Modality:** {c.modality}",
        ])
        if c.sequence:
            lines.append(f"- **Sequence:** `{c.sequence}`")
        if c.smiles:
            lines.append(f"- **SMILES:** `{c.smiles}`")
        if c.iupac_name:
            lines.append(f"- **IUPAC Name:** {c.iupac_name}")
        lines.extend([
            f"- **MW:** {c.molecular_weight:.1f} Da",
            f"- **Net Charge (pH 7.4):** {c.net_charge:+.1f}",
            f"- **GRAVY:** {c.gravy:.3f}",
            f"- **Composite Score:** {c.composite_score:.4f}",
            f"- **Binding Score:** {c.binding_score:.4f}",
            f"- **Selectivity:** {c.selectivity_score:.4f}",
            f"- **Drug-likeness:** {c.drug_likeness:.2f}",
            f"- **ADMET Score:** {c.admet_score:.2f} (flags: {c.admet_flags})",
            f"- **Delivery:** {c.delivery_system}",
            f"- **MOA:** {c.moa_predicted}",
        ])
        if c.sa_score > 0:
            lines.append(f"- **SA Score:** {c.sa_score:.1f}")
        if c.perturbation_score > 0:
            lines.append(f"- **Perturbation Score:** {c.perturbation_score:.3f}")
            if c.cmap_compound_match:
                lines.append(f"- **CMAP Match:** {c.cmap_compound_match}")
        if c.modification:
            lines.append(f"- **Modification:** {c.modification} — {c.modification_detail}")
        lines.append("")

    # Perturbation Biology
    pert_candidates = [c for c in top if c.perturbation_score > 0]
    if pert_candidates:
        lines.extend([
            f"",
            f"## Perturbation Biology Analysis",
            f"",
            f"| ID | CMAP Score | CMAP Match | Network Effect | Disease Reversal | Perturbation |",
            f"|----|-----------|------------|----------------|------------------|-------------|",
        ])
        for c in pert_candidates[:10]:
            lines.append(
                f"| {c.candidate_id} | {c.cmap_connectivity:.3f} | {c.cmap_compound_match or 'N/A'} "
                f"| {c.network_effect_score:.3f} | {c.disease_signature_reversal:.3f} "
                f"| {c.perturbation_score:.3f} |"
            )

    # Synthetic Accessibility
    sa_candidates = [c for c in top if c.sa_score > 0]
    if sa_candidates:
        lines.extend([
            f"",
            f"## Synthetic Accessibility",
            f"",
            f"| ID | SA Score | Interpretation |",
            f"|----|----------|----------------|",
        ])
        for c in sa_candidates[:10]:
            if c.sa_score <= 3:
                interp = "Easy to synthesize"
            elif c.sa_score <= 5:
                interp = "Moderate"
            elif c.sa_score <= 7:
                interp = "Difficult"
            else:
                interp = "Very difficult"
            lines.append(f"| {c.candidate_id} | {c.sa_score:.1f} | {interp} |")

    # Database Source Breakdown
    source_counts: dict[str, int] = {}
    for c in candidates:
        src = c.source or "unknown"
        source_counts[src] = source_counts.get(src, 0) + 1
    if source_counts:
        lines.extend([
            f"",
            f"## Database Source Breakdown",
            f"",
            f"| Source | Count | Percentage |",
            f"|--------|-------|------------|",
        ])
        total = len(candidates) or 1
        for src, count in sorted(source_counts.items(), key=lambda x: -x[1]):
            pct = count / total * 100
            lines.append(f"| {src} | {count} | {pct:.1f}% |")

    # Statistics
    lines.extend([
        f"",
        f"## Pipeline Statistics",
        f"",
        f"| Metric | Value |",
        f"|--------|-------|",
        f"| Total candidates | {len(candidates)} |",
        f"| Library hits | {sum(1 for c in candidates if c.candidate_type == 'library_hit')} |",
        f"| De novo generated | {sum(1 for c in candidates if c.candidate_type == 'de_novo')} |",
        f"| Modified variants | {sum(1 for c in candidates if c.candidate_type == 'modified')} |",
        f"| Database sources | {len(source_counts)} |",
        f"| Mean composite score | {_mean([c.composite_score for c in candidates]):.4f} |",
        f"| Mean ADMET score | {_mean([c.admet_score for c in candidates]):.4f} |",
        f"| Mean SA score | {_mean([c.sa_score for c in sa_candidates]):.1f} |" if sa_candidates else "",
        f"| Candidates with 0 ADMET flags | {sum(1 for c in candidates if c.admet_flags == 0)} |",
        f"",
        f"## Experimental Validation Recommendations",
        f"",
        f"### Tier 1: Immediate Testing (Rank 1-5)",
        f"",
    ])
    for c in top[:5]:
        lines.append(f"- **{c.candidate_id}**: {c.delivery_system} delivery, composite={c.composite_score:.4f}")
    lines.extend([
        f"",
        f"### Tier 2: Secondary Testing (Rank 6-15)",
        f"",
    ])
    for c in top[5:15]:
        lines.append(f"- **{c.candidate_id}**: score={c.composite_score:.4f}")
    lines.extend([
        f"",
        f"### Tier 3: Backup Pool (Rank 16+)",
        f"",
        f"Remaining {max(0, len(candidates) - 15)} candidates available in executive_summary.csv",
        f"",
        f"---",
        f"*Generated by Drug Discovery Pipeline v0.1.0*",
    ])

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# HTML report
# ---------------------------------------------------------------------------

def _candidate_to_json_obj(c: Candidate) -> dict:
    """Export all useful candidate fields for the HTML dashboard."""
    d = c.to_dict()
    d["iupac_name"] = c.iupac_name or ""
    d["gravy"] = round(c.gravy, 3)
    d["modification"] = c.modification or ""
    d["modification_detail"] = c.modification_detail or ""
    d["parent_id"] = c.parent_id or ""
    d["isoelectric_point"] = round(c.isoelectric_point, 2)
    d["structure_confidence"] = round(c.structure_confidence, 3)
    return d


def _build_html_report(
    cfg: PipelineConfig,
    candidates: list[Candidate],
    target: TargetProfile,
) -> str:
    """Build a self-contained interactive HTML dashboard."""
    esc = html.escape
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    gene = esc(target.gene_name)

    # Precompute KPI values
    total = len(candidates)
    n_library = sum(1 for c in candidates if c.candidate_type == 'library_hit')
    n_denovo = sum(1 for c in candidates if c.candidate_type == 'de_novo')
    n_modified = sum(1 for c in candidates if c.candidate_type == 'modified')
    mean_score = _mean([c.composite_score for c in candidates])
    source_set = sorted(set(c.source for c in candidates if c.source))
    n_sources = len(source_set)
    admet_pass = sum(1 for c in candidates if c.admet_flags == 0)
    admet_pct = (admet_pass / total * 100) if total else 0
    has_pert = any(c.perturbation_score > 0 for c in candidates)

    # JSON data blob
    cand_json = json.dumps([_candidate_to_json_obj(c) for c in candidates])
    weights_json = json.dumps(cfg.scoring_weights)
    target_json = json.dumps({
        "gene_name": target.gene_name, "uniprot_id": target.uniprot_id,
        "organism": target.organism, "protein_name": target.protein_name,
        "sequence_length": target.sequence_length, "target_tissue": target.target_tissue,
        "binding_pockets": len(target.binding_pockets),
        "anti_targets": len(target.anti_targets),
    })
    config_json = json.dumps({
        "modality": cfg.modality.value, "mode": cfg.mode.value,
        "device": cfg.device, "top_n": cfg.top_n,
    })

    # Build source option tags
    source_options = "".join(
        f'<option value="{esc(s)}">{esc(s)}</option>' for s in source_set
    )

    # Perturbation tab button (conditional)
    pert_tab_btn = (
        '<button class="tab-btn" data-tab="perturbation">Perturbation Biology</button>'
        if has_pert else ''
    )

    return f'''<!DOCTYPE html>
<html lang="en" data-theme="light">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1.0">
<title>Drug Discovery Dashboard — {gene}</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4"></script>
<style>
:root,[data-theme="light"]{{
  --bg-primary:#f8fafc;--bg-card:#ffffff;--text-primary:#1e293b;--text-secondary:#64748b;
  --accent:#2563eb;--accent-light:#dbeafe;--border:#e2e8f0;--success:#16a34a;
  --warning:#d97706;--danger:#dc2626;--chart-grid:#e2e8f0;--bg-hover:#f1f5f9;
  --bg-detail:#f8fafc;--shadow:0 1px 3px rgba(0,0,0,.1);
}}
[data-theme="dark"]{{
  --bg-primary:#0f172a;--bg-card:#1e293b;--text-primary:#f1f5f9;--text-secondary:#94a3b8;
  --accent:#3b82f6;--accent-light:#1e3a5f;--border:#334155;--success:#22c55e;
  --warning:#f59e0b;--danger:#ef4444;--chart-grid:#334155;--bg-hover:#334155;
  --bg-detail:#0f172a;--shadow:0 1px 3px rgba(0,0,0,.4);
}}
*{{margin:0;padding:0;box-sizing:border-box}}
body{{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;
  background:var(--bg-primary);color:var(--text-primary);line-height:1.6}}
.dashboard{{max-width:1400px;margin:0 auto;padding:1.5rem}}
/* Header */
header{{display:flex;justify-content:space-between;align-items:center;padding:1rem 1.5rem;
  background:var(--bg-card);border-radius:12px;box-shadow:var(--shadow);margin-bottom:1.5rem}}
header h1{{font-size:1.4rem;font-weight:700}}
header h1 span{{color:var(--accent)}}
.header-meta{{color:var(--text-secondary);font-size:.85rem}}
#theme-toggle{{background:var(--accent);color:#fff;border:none;padding:.5rem 1rem;
  border-radius:8px;cursor:pointer;font-size:.85rem;font-weight:500}}
#theme-toggle:hover{{opacity:.85}}
/* KPI cards */
.kpi-grid{{display:grid;grid-template-columns:repeat(4,1fr);gap:1rem;margin-bottom:1.5rem}}
.kpi-card{{background:var(--bg-card);border-radius:12px;padding:1.25rem;box-shadow:var(--shadow);
  border-left:4px solid var(--accent)}}
.kpi-card.green{{border-left-color:var(--success)}}
.kpi-card.orange{{border-left-color:var(--warning)}}
.kpi-card.red{{border-left-color:var(--danger)}}
.kpi-label{{font-size:.8rem;color:var(--text-secondary);text-transform:uppercase;letter-spacing:.5px;font-weight:600}}
.kpi-value{{font-size:1.8rem;font-weight:700;margin:.25rem 0}}
.kpi-sub{{font-size:.75rem;color:var(--text-secondary)}}
/* Tabs */
.tab-bar{{display:flex;gap:.5rem;margin-bottom:1.5rem;overflow-x:auto;padding-bottom:.25rem}}
.tab-btn{{background:var(--bg-card);border:1px solid var(--border);border-radius:8px;
  padding:.6rem 1.2rem;cursor:pointer;font-size:.85rem;font-weight:500;
  color:var(--text-secondary);transition:all .2s;white-space:nowrap}}
.tab-btn:hover{{color:var(--accent);border-color:var(--accent)}}
.tab-btn.active{{background:var(--accent);color:#fff;border-color:var(--accent)}}
.tab-content{{display:none}}
.tab-content.active{{display:block}}
/* Charts grid */
.chart-grid{{display:grid;grid-template-columns:repeat(2,1fr);gap:1.25rem;margin-bottom:1.5rem}}
.chart-card{{background:var(--bg-card);border-radius:12px;padding:1.25rem;box-shadow:var(--shadow)}}
.chart-card h3{{font-size:.95rem;margin-bottom:.75rem;color:var(--text-primary)}}
.chart-card canvas{{width:100%!important;max-height:300px}}
.chart-grid-3{{display:grid;grid-template-columns:repeat(3,1fr);gap:1.25rem;margin-bottom:1.5rem}}
/* Table */
.table-controls{{display:flex;gap:1rem;margin-bottom:1rem;flex-wrap:wrap;align-items:center}}
#search-box{{padding:.5rem .75rem;border:1px solid var(--border);border-radius:8px;
  font-size:.85rem;width:280px;background:var(--bg-card);color:var(--text-primary)}}
#source-filter{{padding:.5rem;border:1px solid var(--border);border-radius:8px;
  font-size:.85rem;background:var(--bg-card);color:var(--text-primary)}}
.data-table-wrap{{overflow-x:auto;background:var(--bg-card);border-radius:12px;box-shadow:var(--shadow)}}
table{{width:100%;border-collapse:collapse;font-size:.82rem}}
th{{background:var(--accent);color:#fff;padding:.6rem .5rem;text-align:left;cursor:pointer;
  white-space:nowrap;user-select:none;position:sticky;top:0}}
th:hover{{opacity:.9}}
th .arrow{{font-size:.7rem;margin-left:.25rem}}
td{{padding:.5rem;border-bottom:1px solid var(--border);white-space:nowrap}}
tr:hover td{{background:var(--bg-hover)}}
tr.expandable{{cursor:pointer}}
.detail-row td{{padding:0;background:var(--bg-detail)}}
.detail-panel{{padding:1rem 1.5rem;display:grid;grid-template-columns:1fr 1fr;gap:.5rem 2rem;
  font-size:.82rem;border-left:3px solid var(--accent)}}
.detail-panel .dp-label{{color:var(--text-secondary);font-weight:600}}
.detail-panel .dp-value{{color:var(--text-primary);word-break:break-all}}
.detail-panel .dp-full{{grid-column:1/-1}}
/* Tier cards */
.tier-grid{{display:grid;grid-template-columns:repeat(3,1fr);gap:1.25rem;margin-top:1rem}}
.tier-card{{background:var(--bg-card);border-radius:12px;padding:1.25rem;box-shadow:var(--shadow)}}
.tier-card h3{{font-size:1rem;margin-bottom:.75rem}}
.tier-card ul{{list-style:none;padding:0}}
.tier-card li{{padding:.35rem 0;border-bottom:1px solid var(--border);font-size:.82rem}}
.tier-card li:last-child{{border:none}}
/* Responsive */
@media(max-width:900px){{
  .kpi-grid{{grid-template-columns:repeat(2,1fr)}}
  .chart-grid,.chart-grid-3,.tier-grid{{grid-template-columns:1fr}}
  .detail-panel{{grid-template-columns:1fr}}
}}
@media(max-width:600px){{
  .kpi-grid{{grid-template-columns:1fr}}
  .dashboard{{padding:.75rem}}
}}
@media print{{
  .tab-bar,#theme-toggle,#search-box,#source-filter{{display:none}}
  .tab-content{{display:block!important}}
  body{{background:#fff;color:#000}}
}}
</style>
</head>
<body>
<div class="dashboard">
<header>
  <div>
    <h1><span>{gene}</span> Drug Discovery Dashboard</h1>
    <div class="header-meta">{esc(target.protein_name)} &middot; {esc(cfg.modality.value)} &middot; {esc(cfg.mode.value)} &middot; {timestamp}</div>
  </div>
  <button id="theme-toggle" onclick="toggleTheme()">Dark Mode</button>
</header>

<div class="kpi-grid">
  <div class="kpi-card">
    <div class="kpi-label">Total Candidates</div>
    <div class="kpi-value">{total}</div>
    <div class="kpi-sub">{n_library} library &middot; {n_denovo} de novo &middot; {n_modified} modified</div>
  </div>
  <div class="kpi-card green">
    <div class="kpi-label">Mean Composite Score</div>
    <div class="kpi-value">{mean_score:.4f}</div>
    <div class="kpi-sub"><div style="background:var(--border);border-radius:4px;height:6px;margin-top:4px">
      <div style="background:var(--success);width:{min(mean_score*100,100):.0f}%;height:100%;border-radius:4px"></div></div></div>
  </div>
  <div class="kpi-card orange">
    <div class="kpi-label">Database Sources</div>
    <div class="kpi-value">{n_sources}</div>
    <div class="kpi-sub">{", ".join(esc(s) for s in source_set[:4])}</div>
  </div>
  <div class="kpi-card" style="border-left-color:{"var(--success)" if admet_pct>=70 else "var(--danger)"}">
    <div class="kpi-label">ADMET Pass Rate</div>
    <div class="kpi-value">{admet_pct:.0f}%</div>
    <div class="kpi-sub">{admet_pass}/{total} with 0 flags</div>
  </div>
</div>

<div class="tab-bar">
  <button class="tab-btn active" data-tab="overview">Overview</button>
  <button class="tab-btn" data-tab="candidates">Candidates</button>
  <button class="tab-btn" data-tab="scoring">Scoring Analysis</button>
  {pert_tab_btn}
  <button class="tab-btn" data-tab="synthesis">Synthesis &amp; Validation</button>
</div>

<!-- TAB 1: Overview -->
<div id="tab-overview" class="tab-content active">
  <div class="chart-grid">
    <div class="chart-card"><h3>Score Distribution</h3><canvas id="chartScoreDist"></canvas></div>
    <div class="chart-card"><h3>Source Breakdown</h3><canvas id="chartSourcePie"></canvas></div>
    <div class="chart-card"><h3>Scoring Weights</h3><canvas id="chartWeightsRadar"></canvas></div>
    <div class="chart-card"><h3>Pipeline Statistics</h3>
      <table style="font-size:.82rem">
        <tr><td>Total candidates</td><td><strong>{total}</strong></td></tr>
        <tr><td>Library hits</td><td>{n_library}</td></tr>
        <tr><td>De novo generated</td><td>{n_denovo}</td></tr>
        <tr><td>Modified variants</td><td>{n_modified}</td></tr>
        <tr><td>Database sources</td><td>{n_sources}</td></tr>
        <tr><td>Mean ADMET score</td><td>{_mean([c.admet_score for c in candidates]):.4f}</td></tr>
        <tr><td>Binding pockets</td><td>{len(target.binding_pockets)}</td></tr>
        <tr><td>Anti-targets</td><td>{len(target.anti_targets)}</td></tr>
      </table>
    </div>
  </div>
</div>

<!-- TAB 2: Candidates -->
<div id="tab-candidates" class="tab-content">
  <div class="table-controls">
    <input type="text" id="search-box" placeholder="Search by ID, name, or SMILES..." oninput="filterTable()">
    <select id="source-filter" onchange="filterTable()">
      <option value="">All Sources</option>
      {source_options}
    </select>
  </div>
  <div class="data-table-wrap">
    <table id="cand-table">
      <thead><tr>
        <th onclick="sortBy('rank')">Rank <span class="arrow"></span></th>
        <th onclick="sortBy('candidate_id')">ID <span class="arrow"></span></th>
        <th onclick="sortBy('iupac_name')">Name <span class="arrow"></span></th>
        <th onclick="sortBy('candidate_type')">Type <span class="arrow"></span></th>
        <th onclick="sortBy('source')">Source <span class="arrow"></span></th>
        <th onclick="sortBy('composite_score')">Score <span class="arrow"></span></th>
        <th onclick="sortBy('binding_score')">Binding <span class="arrow"></span></th>
        <th onclick="sortBy('selectivity_score')">Selectivity <span class="arrow"></span></th>
        <th onclick="sortBy('drug_likeness')">Drug-likeness <span class="arrow"></span></th>
        <th onclick="sortBy('admet_score')">ADMET <span class="arrow"></span></th>
        <th onclick="sortBy('sa_score')">SA <span class="arrow"></span></th>
        <th onclick="sortBy('perturbation_score')">Perturbation <span class="arrow"></span></th>
        <th onclick="sortBy('delivery_system')">Delivery <span class="arrow"></span></th>
        <th onclick="sortBy('moa_predicted')">MOA <span class="arrow"></span></th>
      </tr></thead>
      <tbody id="cand-tbody"></tbody>
    </table>
  </div>
</div>

<!-- TAB 3: Scoring Analysis -->
<div id="tab-scoring" class="tab-content">
  <div class="chart-grid">
    <div class="chart-card"><h3>Selectivity vs Composite Score</h3><canvas id="chartSelVsScore"></canvas></div>
    <div class="chart-card"><h3>ADMET Distribution</h3><canvas id="chartAdmetDist"></canvas></div>
    <div class="chart-card"><h3>Drug-likeness Distribution</h3><canvas id="chartDLDist"></canvas></div>
    <div class="chart-card"><h3>Score Components (Top 10)</h3><canvas id="chartComponents"></canvas></div>
  </div>
</div>

<!-- TAB 4: Perturbation Biology (conditional) -->
<div id="tab-perturbation" class="tab-content">
  <div class="chart-grid">
    <div class="chart-card"><h3>Network Effect vs Disease Reversal</h3><canvas id="chartPertScatter"></canvas></div>
    <div class="chart-card"><h3>Perturbation Candidates</h3>
      <div class="data-table-wrap"><table id="pert-table">
        <thead><tr><th>ID</th><th>CMAP Score</th><th>CMAP Match</th><th>Network Effect</th><th>Disease Reversal</th><th>Perturbation</th></tr></thead>
        <tbody id="pert-tbody"></tbody>
      </table></div>
    </div>
  </div>
</div>

<!-- TAB 5: Synthesis & Validation -->
<div id="tab-synthesis" class="tab-content">
  <div class="chart-grid">
    <div class="chart-card"><h3>SA Score Distribution</h3><canvas id="chartSADist"></canvas></div>
    <div class="chart-card"><h3>SA Interpretation</h3>
      <div class="data-table-wrap"><table>
        <thead><tr><th>ID</th><th>SA Score</th><th>Interpretation</th></tr></thead>
        <tbody id="sa-tbody"></tbody>
      </table></div>
    </div>
  </div>
  <h3 style="margin:1.5rem 0 .75rem;font-size:1.1rem">Experimental Validation Tiers</h3>
  <div class="tier-grid">
    <div class="tier-card" style="border-top:3px solid var(--success)">
      <h3>Tier 1: Immediate Testing</h3><ul id="tier1-list"></ul>
    </div>
    <div class="tier-card" style="border-top:3px solid var(--warning)">
      <h3>Tier 2: Secondary Testing</h3><ul id="tier2-list"></ul>
    </div>
    <div class="tier-card" style="border-top:3px solid var(--text-secondary)">
      <h3>Tier 3: Backup Pool</h3><ul id="tier3-list"></ul>
    </div>
  </div>
</div>

</div><!-- .dashboard -->

<script>
const DATA={{
  candidates:{cand_json},
  target:{target_json},
  config:{config_json},
  weights:{weights_json}
}};
const COLORS=["#2563eb","#16a34a","#d97706","#dc2626","#8b5cf6","#06b6d4","#ec4899","#f97316"];
let sortCol="rank",sortAsc=true,charts={{}};

function isDark(){{return document.documentElement.dataset.theme==="dark"}}
function gridColor(){{return isDark()?"#334155":"#e2e8f0"}}
function textColor(){{return isDark()?"#f1f5f9":"#1e293b"}}
function subColor(){{return isDark()?"#94a3b8":"#64748b"}}

/* Tabs */
document.querySelectorAll(".tab-btn").forEach(btn=>{{
  btn.addEventListener("click",()=>{{
    document.querySelectorAll(".tab-btn").forEach(b=>b.classList.remove("active"));
    document.querySelectorAll(".tab-content").forEach(t=>t.classList.remove("active"));
    btn.classList.add("active");
    document.getElementById("tab-"+btn.dataset.tab).classList.add("active");
  }});
}});

/* Theme toggle */
function toggleTheme(){{
  const el=document.documentElement;
  const dark=el.dataset.theme==="dark";
  el.dataset.theme=dark?"light":"dark";
  document.getElementById("theme-toggle").textContent=dark?"Dark Mode":"Light Mode";
  Object.values(charts).forEach(ch=>{{
    if(!ch)return;
    if(ch.options.scales){{
      Object.values(ch.options.scales).forEach(s=>{{
        if(s.grid)s.grid.color=gridColor();
        if(s.ticks)s.ticks.color=subColor();
        if(s.pointLabels)s.pointLabels.color=subColor();
      }});
    }}
    if(ch.options.plugins&&ch.options.plugins.legend&&ch.options.plugins.legend.labels)
      ch.options.plugins.legend.labels.color=subColor();
    ch.update();
  }});
}}

/* Sort */
function sortBy(col){{
  if(sortCol===col)sortAsc=!sortAsc;else{{sortCol=col;sortAsc=col==="rank"}}
  renderTable();
  document.querySelectorAll("#cand-table th .arrow").forEach((a,i)=>{{a.textContent=""}});
  const ths=document.querySelectorAll("#cand-table th");
  const cols=["rank","candidate_id","iupac_name","candidate_type","source","composite_score",
    "binding_score","selectivity_score","drug_likeness","admet_score","sa_score",
    "perturbation_score","delivery_system","moa_predicted"];
  const idx=cols.indexOf(col);
  if(idx>=0)ths[idx].querySelector(".arrow").textContent=sortAsc?" \\u25B2":" \\u25BC";
}}

/* Filter */
function filterTable(){{
  const q=document.getElementById("search-box").value.toLowerCase();
  const src=document.getElementById("source-filter").value;
  renderTable(q,src);
}}

/* Render candidate table */
function renderTable(query,srcFilter){{
  query=query||"";srcFilter=srcFilter||"";
  let d=[...DATA.candidates];
  if(query)d=d.filter(c=>(c.candidate_id+"|"+c.iupac_name+"|"+c.smiles+"|"+c.sequence).toLowerCase().includes(query));
  if(srcFilter)d=d.filter(c=>c.source===srcFilter);
  d.sort((a,b)=>{{
    let va=a[sortCol],vb=b[sortCol];
    if(typeof va==="string")va=va.toLowerCase();
    if(typeof vb==="string")vb=vb.toLowerCase();
    if(va<vb)return sortAsc?-1:1;
    if(va>vb)return sortAsc?1:-1;
    return 0;
  }});
  const tbody=document.getElementById("cand-tbody");
  tbody.innerHTML="";
  d.forEach(c=>{{
    const name=(c.iupac_name||"").substring(0,25)||(c.sequence||"").substring(0,15)||(c.smiles||"").substring(0,15);
    const tr=document.createElement("tr");
    tr.className="expandable";
    tr.onclick=()=>toggleRow(c.candidate_id,tr);
    tr.innerHTML=`<td>${{c.rank}}</td><td>${{esc(c.candidate_id)}}</td><td title="${{esc(c.iupac_name||"")}}">${{esc(name)}}</td>`+
      `<td>${{c.candidate_type}}</td><td>${{esc(c.source)}}</td>`+
      `<td><strong>${{c.composite_score.toFixed(4)}}</strong></td>`+
      `<td>${{c.binding_score.toFixed(4)}}</td><td>${{c.selectivity_score.toFixed(4)}}</td>`+
      `<td>${{c.drug_likeness.toFixed(2)}}</td><td>${{c.admet_score.toFixed(2)}}</td>`+
      `<td>${{c.sa_score>0?c.sa_score.toFixed(1):"-"}}</td>`+
      `<td>${{c.perturbation_score>0?c.perturbation_score.toFixed(3):"-"}}</td>`+
      `<td>${{esc(c.delivery_system)}}</td><td>${{esc(c.moa_predicted)}}</td>`;
    tbody.appendChild(tr);
  }});
}}
function esc(s){{const d=document.createElement("div");d.textContent=s||"";return d.innerHTML}}

/* Toggle detail row */
function toggleRow(id,tr){{
  const next=tr.nextElementSibling;
  if(next&&next.classList.contains("detail-row")){{next.remove();return}}
  const c=DATA.candidates.find(x=>x.candidate_id===id);
  if(!c)return;
  const dr=document.createElement("tr");
  dr.className="detail-row";
  const smiles=c.smiles||c.sequence||"N/A";
  let extra="";
  if(c.modification)extra+=`<div class="dp-label">Modification</div><div class="dp-value">${{esc(c.modification)}} — ${{esc(c.modification_detail)}}</div>`;
  if(c.parent_id)extra+=`<div class="dp-label">Parent ID</div><div class="dp-value">${{esc(c.parent_id)}}</div>`;
  if(c.cmap_compound_match)extra+=`<div class="dp-label">CMAP Match</div><div class="dp-value">${{esc(c.cmap_compound_match)}}</div>`;
  dr.innerHTML=`<td colspan="14"><div class="detail-panel">
    <div class="dp-label">SMILES/Sequence</div><div class="dp-value dp-full" style="grid-column:1/-1;font-family:monospace;font-size:.75rem">${{esc(smiles)}}</div>
    <div class="dp-label">Molecular Weight</div><div class="dp-value">${{c.molecular_weight.toFixed(1)}} Da</div>
    <div class="dp-label">Net Charge (pH 7.4)</div><div class="dp-value">${{c.net_charge>=0?"+":"" }}${{c.net_charge.toFixed(1)}}</div>
    <div class="dp-label">GRAVY</div><div class="dp-value">${{c.gravy.toFixed(3)}}</div>
    <div class="dp-label">Isoelectric Point</div><div class="dp-value">${{c.isoelectric_point.toFixed(2)}}</div>
    <div class="dp-label">Binding Score</div><div class="dp-value">${{c.binding_score.toFixed(4)}}</div>
    <div class="dp-label">Selectivity</div><div class="dp-value">${{c.selectivity_score.toFixed(4)}}</div>
    <div class="dp-label">Drug-likeness</div><div class="dp-value">${{c.drug_likeness.toFixed(4)}}</div>
    <div class="dp-label">ADMET Score</div><div class="dp-value">${{c.admet_score.toFixed(4)}} (${{c.admet_flags}} flags)</div>
    <div class="dp-label">SA Score</div><div class="dp-value">${{c.sa_score>0?c.sa_score.toFixed(1):"N/A"}}</div>
    <div class="dp-label">Perturbation</div><div class="dp-value">${{c.perturbation_score>0?c.perturbation_score.toFixed(3):"N/A"}}</div>
    <div class="dp-label">CMAP Connectivity</div><div class="dp-value">${{c.cmap_connectivity.toFixed(3)}}</div>
    <div class="dp-label">Network Effect</div><div class="dp-value">${{c.network_effect_score.toFixed(3)}}</div>
    <div class="dp-label">Disease Reversal</div><div class="dp-value">${{c.disease_signature_reversal.toFixed(3)}}</div>
    <div class="dp-label">Delivery</div><div class="dp-value">${{esc(c.delivery_system)}}</div>
    <div class="dp-label">MOA</div><div class="dp-value">${{esc(c.moa_predicted)}}</div>
    ${{extra}}
  </div></td>`;
  tr.after(dr);
}}

/* Histogram helper */
function makeHistData(values,bins,color){{
  if(!values.length)return{{labels:[],data:[]}};
  const mn=Math.min(...values),mx=Math.max(...values);
  const step=(mx-mn)/bins||1;
  const labels=[],counts=new Array(bins).fill(0);
  for(let i=0;i<bins;i++)labels.push((mn+step*i).toFixed(2));
  values.forEach(v=>{{let idx=Math.min(Math.floor((v-mn)/step),bins-1);if(idx<0)idx=0;counts[idx]++}});
  return{{labels,data:counts}};
}}

/* Init charts */
function initCharts(){{
  const cands=DATA.candidates;
  // 1. Score distribution
  const scores=cands.map(c=>c.composite_score);
  const sd=makeHistData(scores,15);
  charts.scoreDist=new Chart(document.getElementById("chartScoreDist"),{{
    type:"bar",data:{{labels:sd.labels,datasets:[{{label:"Count",data:sd.data,backgroundColor:"#2563eb"}}]}},
    options:{{responsive:true,plugins:{{legend:{{display:false}}}},scales:{{x:{{title:{{display:true,text:"Composite Score",color:subColor()}},grid:{{color:gridColor()}},ticks:{{color:subColor()}}}},y:{{title:{{display:true,text:"Count",color:subColor()}},grid:{{color:gridColor()}},ticks:{{color:subColor()}}}}}}}}
  }});

  // 2. Source breakdown
  const srcMap={{}};cands.forEach(c=>{{srcMap[c.source]=(srcMap[c.source]||0)+1}});
  const srcLabels=Object.keys(srcMap),srcVals=Object.values(srcMap);
  charts.sourcePie=new Chart(document.getElementById("chartSourcePie"),{{
    type:"doughnut",data:{{labels:srcLabels,datasets:[{{data:srcVals,backgroundColor:COLORS.slice(0,srcLabels.length)}}]}},
    options:{{responsive:true,plugins:{{legend:{{position:"right",labels:{{color:subColor()}}}}}}}}
  }});

  // 3. Weights radar
  const wLabels=Object.keys(DATA.weights),wVals=Object.values(DATA.weights);
  charts.weightsRadar=new Chart(document.getElementById("chartWeightsRadar"),{{
    type:"radar",data:{{labels:wLabels,datasets:[{{label:"Weight",data:wVals,backgroundColor:"rgba(37,99,235,.2)",borderColor:"#2563eb",pointBackgroundColor:"#2563eb"}}]}},
    options:{{responsive:true,scales:{{r:{{beginAtZero:true,grid:{{color:gridColor()}},pointLabels:{{color:subColor()}},ticks:{{color:subColor(),backdropColor:"transparent"}}}}}},plugins:{{legend:{{labels:{{color:subColor()}}}}}}}}
  }});

  // 4. Selectivity vs Score scatter
  const srcList=[...new Set(cands.map(c=>c.source))];
  const scatterDS=srcList.map((src,i)=>{{
    const pts=cands.filter(c=>c.source===src).map(c=>({{x:c.selectivity_score,y:c.composite_score}}));
    return{{label:src,data:pts,backgroundColor:COLORS[i%COLORS.length],pointRadius:5}};
  }});
  charts.selVsScore=new Chart(document.getElementById("chartSelVsScore"),{{
    type:"scatter",data:{{datasets:scatterDS}},
    options:{{responsive:true,scales:{{x:{{title:{{display:true,text:"Selectivity",color:subColor()}},grid:{{color:gridColor()}},ticks:{{color:subColor()}}}},y:{{title:{{display:true,text:"Composite Score",color:subColor()}},grid:{{color:gridColor()}},ticks:{{color:subColor()}}}}}},plugins:{{legend:{{labels:{{color:subColor()}}}}}}}}
  }});

  // 5. ADMET distribution
  const ad=makeHistData(cands.map(c=>c.admet_score),12);
  charts.admetDist=new Chart(document.getElementById("chartAdmetDist"),{{
    type:"bar",data:{{labels:ad.labels,datasets:[{{label:"Count",data:ad.data,backgroundColor:"#d97706"}}]}},
    options:{{responsive:true,plugins:{{legend:{{display:false}}}},scales:{{x:{{title:{{display:true,text:"ADMET Score",color:subColor()}},grid:{{color:gridColor()}},ticks:{{color:subColor()}}}},y:{{grid:{{color:gridColor()}},ticks:{{color:subColor()}}}}}}}}
  }});

  // 6. Drug-likeness distribution
  const dl=makeHistData(cands.map(c=>c.drug_likeness),12);
  charts.dlDist=new Chart(document.getElementById("chartDLDist"),{{
    type:"bar",data:{{labels:dl.labels,datasets:[{{label:"Count",data:dl.data,backgroundColor:"#16a34a"}}]}},
    options:{{responsive:true,plugins:{{legend:{{display:false}}}},scales:{{x:{{title:{{display:true,text:"Drug-likeness",color:subColor()}},grid:{{color:gridColor()}},ticks:{{color:subColor()}}}},y:{{grid:{{color:gridColor()}},ticks:{{color:subColor()}}}}}}}}
  }});

  // 7. Score components stacked bar (top 10)
  const top10=cands.slice(0,10);
  const compLabels=top10.map(c=>c.candidate_id);
  const compKeys=["binding_score","selectivity_score","drug_likeness","admet_score","perturbation_score"];
  const compColors=["#2563eb","#8b5cf6","#16a34a","#d97706","#ec4899"];
  const compDS=compKeys.map((k,i)=>({{
    label:k.replace("_"," "),data:top10.map(c=>c[k]||0),backgroundColor:compColors[i]
  }}));
  charts.components=new Chart(document.getElementById("chartComponents"),{{
    type:"bar",data:{{labels:compLabels,datasets:compDS}},
    options:{{responsive:true,indexAxis:"y",scales:{{x:{{stacked:true,grid:{{color:gridColor()}},ticks:{{color:subColor()}}}},y:{{stacked:true,grid:{{color:gridColor()}},ticks:{{color:subColor(),font:{{size:10}}}}}}}},plugins:{{legend:{{labels:{{color:subColor(),font:{{size:10}}}}}}}}}}
  }});

  // 8. SA Score distribution (if any)
  const saVals=cands.map(c=>c.sa_score).filter(v=>v>0);
  if(saVals.length){{
    const sa=makeHistData(saVals,10);
    charts.saDist=new Chart(document.getElementById("chartSADist"),{{
      type:"bar",data:{{labels:sa.labels,datasets:[{{label:"Count",data:sa.data,backgroundColor:"#8b5cf6"}}]}},
      options:{{responsive:true,plugins:{{legend:{{display:false}}}},scales:{{x:{{title:{{display:true,text:"SA Score (1=easy, 10=hard)",color:subColor()}},grid:{{color:gridColor()}},ticks:{{color:subColor()}}}},y:{{grid:{{color:gridColor()}},ticks:{{color:subColor()}}}}}}}}
    }});
  }}

  // Perturbation scatter
  const pertCands=cands.filter(c=>c.perturbation_score>0);
  if(pertCands.length){{
    charts.pertScatter=new Chart(document.getElementById("chartPertScatter"),{{
      type:"scatter",data:{{datasets:[{{
        label:"Candidates",
        data:pertCands.map(c=>({{x:c.network_effect_score,y:c.disease_signature_reversal}})),
        backgroundColor:"#ec4899",pointRadius:6
      }}]}},
      options:{{responsive:true,scales:{{x:{{title:{{display:true,text:"Network Effect",color:subColor()}},grid:{{color:gridColor()}},ticks:{{color:subColor()}}}},y:{{title:{{display:true,text:"Disease Reversal",color:subColor()}},grid:{{color:gridColor()}},ticks:{{color:subColor()}}}}}},plugins:{{legend:{{display:false}}}}}}
    }});
    // Perturbation table
    const ptb=document.getElementById("pert-tbody");
    pertCands.forEach(c=>{{
      ptb.innerHTML+=`<tr><td>${{esc(c.candidate_id)}}</td><td>${{c.cmap_connectivity.toFixed(3)}}</td><td>${{esc(c.cmap_compound_match||"N/A")}}</td><td>${{c.network_effect_score.toFixed(3)}}</td><td>${{c.disease_signature_reversal.toFixed(3)}}</td><td>${{c.perturbation_score.toFixed(3)}}</td></tr>`;
    }});
  }}

  // SA table
  const saCands=cands.filter(c=>c.sa_score>0);
  const satb=document.getElementById("sa-tbody");
  saCands.forEach(c=>{{
    let interp=c.sa_score<=3?"Easy":c.sa_score<=5?"Moderate":c.sa_score<=7?"Difficult":"Very difficult";
    satb.innerHTML+=`<tr><td>${{esc(c.candidate_id)}}</td><td>${{c.sa_score.toFixed(1)}}</td><td>${{interp}}</td></tr>`;
  }});

  // Tier lists
  const sorted=[...cands].sort((a,b)=>a.rank-b.rank);
  const t1=document.getElementById("tier1-list"),t2=document.getElementById("tier2-list"),t3=document.getElementById("tier3-list");
  sorted.slice(0,5).forEach(c=>{{t1.innerHTML+=`<li><strong>${{esc(c.candidate_id)}}</strong> — ${{esc(c.delivery_system)}}, score=${{c.composite_score.toFixed(4)}}</li>`}});
  sorted.slice(5,15).forEach(c=>{{t2.innerHTML+=`<li><strong>${{esc(c.candidate_id)}}</strong> — score=${{c.composite_score.toFixed(4)}}</li>`}});
  const remaining=Math.max(0,sorted.length-15);
  t3.innerHTML=`<li>${{remaining}} candidates in executive_summary.csv</li>`;
}}

document.addEventListener("DOMContentLoaded",()=>{{
  renderTable();
  initCharts();
}});
</script>
</body>
</html>'''


# ---------------------------------------------------------------------------
# Visualizations
# ---------------------------------------------------------------------------

def _generate_figures(candidates: list[Candidate], fig_dir: Path) -> None:
    """Generate analysis plots."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        logger.warning("matplotlib/seaborn not available; skipping figures")
        return

    scores = [c.composite_score for c in candidates]
    dl_scores = [c.drug_likeness for c in candidates]
    sel_scores = [c.selectivity_score for c in candidates]
    admet_scores = [c.admet_score for c in candidates]

    # Composite score distribution
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(scores, bins=20, color="#2980b9", edgecolor="white", alpha=0.8)
    ax.set_xlabel("Composite Score")
    ax.set_ylabel("Count")
    ax.set_title("Composite Score Distribution")
    fig.tight_layout()
    fig.savefig(fig_dir / "composite_score_distribution.png", dpi=150)
    plt.close(fig)

    # Drug-likeness distribution
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(dl_scores, bins=20, color="#27ae60", edgecolor="white", alpha=0.8)
    ax.set_xlabel("Drug-likeness Score")
    ax.set_ylabel("Count")
    ax.set_title("Drug-likeness Distribution")
    fig.tight_layout()
    fig.savefig(fig_dir / "drug_likeness_distribution.png", dpi=150)
    plt.close(fig)

    # Selectivity vs composite score scatter
    fig, ax = plt.subplots(figsize=(8, 5))
    sources = list(set(c.source for c in candidates))
    colors = plt.cm.Set2.colors
    for i, src in enumerate(sources):
        src_cands = [c for c in candidates if c.source == src]
        ax.scatter(
            [c.selectivity_score for c in src_cands],
            [c.composite_score for c in src_cands],
            label=src, color=colors[i % len(colors)], alpha=0.7, s=40,
        )
    ax.set_xlabel("Selectivity Score")
    ax.set_ylabel("Composite Score")
    ax.set_title("Selectivity vs Composite Score")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(fig_dir / "selectivity_vs_score.png", dpi=150)
    plt.close(fig)

    # ADMET score distribution
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(admet_scores, bins=20, color="#e67e22", edgecolor="white", alpha=0.8)
    ax.set_xlabel("ADMET Score")
    ax.set_ylabel("Count")
    ax.set_title("ADMET Score Distribution")
    fig.tight_layout()
    fig.savefig(fig_dir / "admet_score_distribution.png", dpi=150)
    plt.close(fig)

    # Source comparison (box plot)
    if len(sources) > 1:
        fig, ax = plt.subplots(figsize=(8, 5))
        data_by_source = {src: [c.composite_score for c in candidates if c.source == src] for src in sources}
        ax.boxplot(data_by_source.values(), labels=data_by_source.keys())
        ax.set_ylabel("Composite Score")
        ax.set_title("Score by Source")
        fig.tight_layout()
        fig.savefig(fig_dir / "source_comparison.png", dpi=150)
        plt.close(fig)

    # Database source breakdown (pie chart)
    if len(sources) > 1:
        fig, ax = plt.subplots(figsize=(7, 7))
        source_counts = {src: sum(1 for c in candidates if c.source == src) for src in sources}
        labels = list(source_counts.keys())
        sizes = list(source_counts.values())
        ax.pie(sizes, labels=labels, autopct="%1.1f%%", startangle=140)
        ax.set_title("Candidates by Database Source")
        fig.tight_layout()
        fig.savefig(fig_dir / "source_breakdown.png", dpi=150)
        plt.close(fig)

    # SA score distribution (small molecules only)
    sa_scores = [c.sa_score for c in candidates if c.sa_score > 0]
    if sa_scores:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.hist(sa_scores, bins=15, color="#8e44ad", edgecolor="white", alpha=0.8)
        ax.set_xlabel("SA Score (1=easy, 10=hard)")
        ax.set_ylabel("Count")
        ax.set_title("Synthetic Accessibility Distribution")
        fig.tight_layout()
        fig.savefig(fig_dir / "sa_score_distribution.png", dpi=150)
        plt.close(fig)

    logger.info("Figures saved to %s", fig_dir)


def _mean(values: list[float]) -> float:
    """Safe mean calculation."""
    if not values:
        return 0.0
    return sum(values) / len(values)
