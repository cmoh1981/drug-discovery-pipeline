"""M9: Final report generation.

Produces HTML, Markdown, and CSV reports with visualizations.
Refactored from YARS2 pipeline: analyze_results.py
"""

from __future__ import annotations

import html
import logging
from datetime import datetime
from pathlib import Path
from typing import Optional

from drugdiscovery.types import Candidate, PipelineConfig, TargetProfile
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
        header = "| Rank | ID | Type | Source | Score | Drug-likeness | ADMET | Delivery |"
        sep = "|------|-----|------|--------|-------|---------------|-------|----------|"
        lines.append(header)
        lines.append(sep)
        for c in top:
            seq_or_smiles = c.sequence[:20] if c.sequence else (c.smiles[:20] if c.smiles else "")
            lines.append(
                f"| {c.rank} | {c.candidate_id} | {c.candidate_type} | {c.source} "
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
        if c.modification:
            lines.append(f"- **Modification:** {c.modification} — {c.modification_detail}")
        lines.append("")

    # Statistics
    lines.extend([
        f"## Pipeline Statistics",
        f"",
        f"| Metric | Value |",
        f"|--------|-------|",
        f"| Total candidates | {len(candidates)} |",
        f"| Library hits | {sum(1 for c in candidates if c.candidate_type == 'library_hit')} |",
        f"| De novo generated | {sum(1 for c in candidates if c.candidate_type == 'de_novo')} |",
        f"| Modified variants | {sum(1 for c in candidates if c.candidate_type == 'modified')} |",
        f"| Mean composite score | {_mean([c.composite_score for c in candidates]):.4f} |",
        f"| Mean ADMET score | {_mean([c.admet_score for c in candidates]):.4f} |",
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

def _build_html_report(
    cfg: PipelineConfig,
    candidates: list[Candidate],
    target: TargetProfile,
) -> str:
    """Build a self-contained HTML report."""
    top = candidates[:cfg.top_n]
    rows_html = ""
    for c in top:
        seq_display = c.sequence[:15] + '...' if len(c.sequence) > 15 else c.sequence or c.smiles[:15]
        rows_html += f"""<tr>
            <td>{c.rank}</td>
            <td>{html.escape(str(c.candidate_id))}</td>
            <td>{html.escape(str(c.source))}</td>
            <td>{html.escape(seq_display)}</td>
            <td>{c.composite_score:.4f}</td>
            <td>{c.drug_likeness:.2f}</td>
            <td>{c.admet_score:.2f}</td>
            <td>{html.escape(str(c.delivery_system))}</td>
        </tr>"""

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Drug Discovery Report - {html.escape(target.gene_name)}</title>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; margin: 2em; color: #333; max-width: 1200px; margin: 0 auto; padding: 2em; }}
h1 {{ color: #1a5276; border-bottom: 3px solid #2980b9; padding-bottom: 0.5em; }}
h2 {{ color: #2c3e50; margin-top: 2em; }}
table {{ border-collapse: collapse; width: 100%; margin: 1em 0; }}
th, td {{ border: 1px solid #ddd; padding: 8px 12px; text-align: left; }}
th {{ background-color: #2980b9; color: white; }}
tr:nth-child(even) {{ background-color: #f8f9fa; }}
tr:hover {{ background-color: #e8f4fd; }}
.metric {{ display: inline-block; background: #ecf0f1; padding: 0.5em 1em; margin: 0.3em; border-radius: 4px; }}
.metric strong {{ color: #2980b9; }}
.summary {{ background: #eaf2f8; padding: 1em; border-radius: 8px; margin: 1em 0; }}
img {{ max-width: 100%; height: auto; margin: 1em 0; }}
</style>
</head>
<body>
<h1>Drug Discovery Pipeline Report</h1>
<p><em>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</em></p>

<div class="summary">
<h2>Target: {html.escape(target.gene_name)} ({html.escape(target.uniprot_id)})</h2>
<div class="metric"><strong>Organism:</strong> {html.escape(target.organism)}</div>
<div class="metric"><strong>Length:</strong> {target.sequence_length} aa</div>
<div class="metric"><strong>Modality:</strong> {cfg.modality.value}</div>
<div class="metric"><strong>Mode:</strong> {cfg.mode.value}</div>
<div class="metric"><strong>Candidates:</strong> {len(candidates)}</div>
</div>

<h2>Top {len(top)} Candidates</h2>
<table>
<thead>
<tr><th>Rank</th><th>ID</th><th>Source</th><th>Sequence/SMILES</th><th>Score</th><th>Drug-likeness</th><th>ADMET</th><th>Delivery</th></tr>
</thead>
<tbody>
{rows_html}
</tbody>
</table>

<h2>Score Distribution</h2>
<p><img src="figures/composite_score_distribution.png" alt="Score Distribution"></p>

<h2>Figures</h2>
<p><img src="figures/drug_likeness_distribution.png" alt="Drug-likeness Distribution"></p>
<p><img src="figures/selectivity_vs_score.png" alt="Selectivity vs Score"></p>

<hr>
<p><em>Drug Discovery Pipeline v0.1.0</em></p>
</body>
</html>"""
    return html


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

    logger.info("Figures saved to %s", fig_dir)


def _mean(values: list[float]) -> float:
    """Safe mean calculation."""
    if not values:
        return 0.0
    return sum(values) / len(values)
