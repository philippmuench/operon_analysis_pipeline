#!/usr/bin/env python3
"""Visualize codon-level variation for operon alignments.

This helper script reproduces the simple codon counting logic used in the
selection analysis and emits per-codon summaries plus an annotated plot that
highlights which positions are variable and how they vary (synonymous vs
non-synonymous). The goal is to provide ready-to-plot material for figures.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib as mpl
from collections import Counter

import matplotlib.pyplot as plt
import pandas as pd
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import numpy as np

CURRENT_DIR = Path(__file__).resolve().parent
ROOT_DIR = CURRENT_DIR.parent
CODON_PLOT_DIR = ROOT_DIR / '14_codon_plots'
if str(CODON_PLOT_DIR) not in sys.path:
    sys.path.append(str(CODON_PLOT_DIR))

from codon_window_plots import (
    compute_codon_stats,
    plot_codon_window_codons,
    plot_codon_window_combined,
)

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
# Standard bacterial genetic code (NCBI table 11) without stop codons.
GENETIC_CODE: Dict[str, str] = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}
GENETIC_CODE = {codon: aa for codon, aa in GENETIC_CODE.items() if aa != '*'}

CATEGORY_COLORS = {
    'synonymous-only': '#2E7D32',
    'nonsynonymous-only': '#C62828',
}

NUCLEOTIDE_COLORS = {
    'A': '#8dd3c7',
    'C': '#ffffb3',
    'G': '#bebada',
    'T': '#fb8072'
}


@dataclass
class CodonSummary:
    """Per-codon summary with the same logic as the analysis pipeline."""

    codon_index: int
    n_sequences: int
    n_variants: int
    syn_changes: int
    nonsyn_changes: int
    category: str
    aa_states: Dict[str, int]
    codon_counts: Dict[str, int]

    @property
    def most_common_codons(self) -> str:
        items = sorted(self.codon_counts.items(), key=lambda kv: (-kv[1], kv[0]))
        return ", ".join(f"{codon}×{count}" for codon, count in items)

    @property
    def aa_state_string(self) -> str:
        items = sorted(self.aa_states.items(), key=lambda kv: (-kv[1], kv[0]))
        return ", ".join(f"{aa}×{count}" for aa, count in items)


def iter_variable_codons(alignment) -> Iterable[CodonSummary]:
    seq_length = alignment.get_alignment_length()
    n_codons = seq_length // 3

    for codon_idx in range(n_codons):
        start = codon_idx * 3
        codons: List[str] = []

        for record in alignment:
            codon = str(record.seq[start:start + 3]).upper()
            if len(codon) != 3 or '-' in codon or 'N' in codon:
                continue
            if codon not in GENETIC_CODE:
                continue
            codons.append(codon)

        if len(codons) < 2:
            continue

        codon_counts = Counter(codons)
        if len(codon_counts) < 2:
            continue

        aa_counts = Counter(GENETIC_CODE[c] for c in codon_counts)
        syn = int(len(aa_counts) == 1)
        nonsyn = int(len(aa_counts) > 1)

        if syn:
            category = 'synonymous-only'
        else:
            category = 'nonsynonymous-only'

        yield CodonSummary(
            codon_index=codon_idx + 1,
            n_sequences=sum(codon_counts.values()),
            n_variants=len(codon_counts),
            syn_changes=syn,
            nonsyn_changes=nonsyn,
            category=category,
            aa_states=dict(aa_counts),
            codon_counts=dict(codon_counts),
        )


def summarise_gene(alignment_path: Path) -> Tuple[pd.DataFrame, List[CodonSummary], MultipleSeqAlignment]:
    alignment = AlignIO.read(alignment_path, 'fasta')
    records = list(iter_variable_codons(alignment))
    if not records:
        return pd.DataFrame(), [], alignment

    data = [
        {
            'codon_index': rec.codon_index,
            'category': rec.category,
            'syn_changes': rec.syn_changes,
            'nonsyn_changes': rec.nonsyn_changes,
            'n_variants': rec.n_variants,
            'n_sequences': rec.n_sequences,
            'aa_states': rec.aa_state_string,
            'codon_counts': rec.most_common_codons,
        }
        for rec in records
    ]

    df = pd.DataFrame(data).sort_values('codon_index').reset_index(drop=True)
    return df, records, alignment


def plot_gene(df: pd.DataFrame, gene: str, output_dir: Path,
              highlight_ranges: list[tuple[float, float]] | None = None) -> None:
    if df.empty:
        return

    df = df.reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(12, 4.5))

    df = df.reset_index(drop=True)

    diversity = df['aa_states'].apply(lambda s: len(str(s).split(', ')) if pd.notna(s) else 0)
    y = diversity
    colors = [CATEGORY_COLORS.get(cat, '#555555') for cat in df['category']]
    sizes = [60 + 40 * max(0, n - 2) for n in df['n_variants']]  # Highlight diverse codons

    if highlight_ranges:
        for start, end in highlight_ranges:
            ax.axvspan(start, end, color='#ffeb3b', alpha=0.2, zorder=0)

    ax.scatter(df['codon_index'], y, c=colors, s=sizes, linewidths=0.6,
               edgecolors='black', alpha=0.55, zorder=1)

    # Label sites where more than two unique codons occur
    for idx in df[df['n_variants'] > 2].index:
        row = df.loc[idx]
        ax.text(row['codon_index'], y[idx], str(int(row['n_variants'])),
                fontsize=8, fontweight='bold', color='white', ha='center', va='center')

    ax.annotate('Numbers indicate distinct codon alleles at that position',
                xy=(0.02, 0.98), xycoords='axes fraction', fontsize=9,
                ha='left', va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.6))

    ax.set_xlabel('Codon position')
    max_y = max(y) if len(y) else 0
    ax.set_ylabel('Distinct amino acids')
    ax.set_ylim(-0.5, max_y + 1.5)
    min_x = df['codon_index'].min() - 1
    max_x = df['codon_index'].max() + 1
    if highlight_ranges:
        min_x = min(min_x, min(r[0] for r in highlight_ranges) - 0.5)
        max_x = max(max_x, max(r[1] for r in highlight_ranges) + 0.5)
    ax.set_xlim(min_x, max_x)
    ax.grid(axis='x', alpha=0.2)

    # Legend outside plot area
    legend_handles = [
        plt.Line2D([], [], marker='o', linestyle='', color=CATEGORY_COLORS.get(key, '#555555'),
                   markeredgecolor='black', markersize=8, label=label)
        for key, label in [('synonymous-only', 'Synonymous'), ('nonsynonymous-only', 'Non-synonymous')]
    ]
    ax.legend(handles=legend_handles, loc='upper left', bbox_to_anchor=(1.02, 1),
              borderaxespad=0, frameon=False, title='Category')

    fig.tight_layout()
    out_path = output_dir / f'{gene}_variable_codons.pdf'
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists():
        out_path.unlink()
    fig.savefig(out_path)
    plt.close(fig)


def plot_focus_codons(records: List[CodonSummary], alignment: MultipleSeqAlignment,
                      focus_codons: Iterable[int], gene: str, output_dir: Path,
                      max_examples: int = 3, window_size: int = 5) -> None:
    if not focus_codons:
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    example_dir = output_dir / 'examples'
    example_dir.mkdir(parents=True, exist_ok=True)
    window_dir = output_dir / 'windows'
    window_dir.mkdir(parents=True, exist_ok=True)
    palette = plt.get_cmap('tab20')
    total_codons = alignment.get_alignment_length() // 3

    for codon_idx in focus_codons:
        idx = int(codon_idx)
        if idx < 1 or idx > total_codons:
            print(f'[warning] Focus codon {idx} out of range for {gene} (1-{total_codons})')
            continue

        stats = compute_codon_stats(alignment, idx)
        codon_counts = stats['codon_counts']
        if not codon_counts:
            print(f'[warning] Focus codon {idx} has no high-quality bases for {gene}')
            continue

        if stats['category'] == 'invariant':
            print(f'[info] Focus codon {idx} is invariant for {gene}; generating context plots.')

        items = sorted(codon_counts.items(), key=lambda kv: (-kv[1], kv[0]))
        aa_counts = stats['aa_counts']
        unique_aas = sorted(aa_counts, key=lambda aa: (-aa_counts[aa], aa))
        aa_colors = {aa: palette(i % palette.N) for i, aa in enumerate(unique_aas)}

        plot_alignment_examples(alignment, idx, items, aa_colors, gene,
                                 example_dir, max_examples)
        plot_codon_window_heatmap(alignment, idx, gene, window_dir,
                                  window_size)


def plot_amino_acid_diversity(records: List[CodonSummary], gene: str,
                              output_dir: Path) -> None:
    if not records:
        return

    counts = Counter(len(rec.aa_states) for rec in records if rec.aa_states)
    if not counts:
        return

    xs = sorted(counts)
    values = [counts[x] for x in xs]

    fig, ax = plt.subplots(figsize=(3, 4))
    bar_width = 0.4
    bars = ax.bar([x - bar_width/2 for x in xs], values, width=bar_width,
                  color='#377eb8', edgecolor='black', alpha=0.8)
    ax.set_xlabel('Distinct amino acids per codon')
    ax.set_ylabel('Number of codons')
    ax.set_title(f'{gene}: amino-acid diversity in variable codons')
    ax.set_xticks(xs)
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                str(value), ha='center', va='bottom', fontsize=9, fontweight='bold')

    fig.tight_layout()
    out_path = output_dir / f'{gene}_aa_diversity.pdf'
    if out_path.exists():
        out_path.unlink()
    fig.savefig(out_path)
    plt.close(fig)


def plot_alignment_examples(alignment: MultipleSeqAlignment, codon_idx: int,
                            codon_items: List[Tuple[str, int]],
                            aa_colors: Dict[str, str], gene: str,
                            output_dir: Path, max_examples: int) -> None:
    start = (codon_idx - 1) * 3
    end = start + 3

    if start < 0 or end > alignment.get_alignment_length():
        return

    codon_to_sequences: Dict[str, List[str]] = {}
    for record in alignment:
        codon = str(record.seq[start:end]).upper()
        if len(codon) != 3 or '-' in codon or 'N' in codon or codon not in GENETIC_CODE:
            continue
        codon_to_sequences.setdefault(codon, []).append(record.id)

    rows: List[Tuple[str, str, str]] = []
    for codon, _ in codon_items:
        seq_ids = codon_to_sequences.get(codon, [])[:max_examples]
        for seq_id in seq_ids:
            rows.append((seq_id, codon, GENETIC_CODE[codon]))

    if not rows:
        return

    cell_text = []
    cell_colors = []
    row_labels = []
    for seq_id, codon, aa in rows:
        row_labels.append(seq_id[:25])
        nucleotides = list(codon)
        cell_text.append(nucleotides + [aa])
        cell_colors.append(
            [NUCLEOTIDE_COLORS.get(n, '#e0e0e0') for n in nucleotides] +
            [aa_colors.get(aa, '#d9d9d9')]
        )

    fig_height = max(2.0, 1 + 0.4 * len(rows))
    fig, ax = plt.subplots(figsize=(5, fig_height))
    ax.axis('off')

    table = ax.table(
        cellText=cell_text,
        cellColours=cell_colors,
        rowLabels=row_labels,
        colLabels=['base1', 'base2', 'base3', 'AA'],
        cellLoc='center',
        loc='center'
    )
    table.scale(1, 1.2)
    table.auto_set_font_size(False)
    table.set_fontsize(9)

    ax.set_title(
        f'{gene} codon {codon_idx} allele examples',
        fontweight='bold', fontsize=12
    )

    out_path = output_dir / f'{gene}_codon_{codon_idx}_examples.pdf'
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
def plot_codon_window_heatmap(alignment: MultipleSeqAlignment, focus_idx: int,
                              gene: str, output_dir: Path, window_size: int) -> None:
    total_codons = alignment.get_alignment_length() // 3
    window_min = max(1, focus_idx - window_size)
    window_max = min(total_codons, focus_idx + window_size)
    positions = list(range(window_min, window_max + 1))

    stats_per_pos = {pos: compute_codon_stats(alignment, pos) for pos in positions}
    focus_col = positions.index(focus_idx)

    amino_acids = []
    for stats in stats_per_pos.values():
        for aa in stats['aa_counts']:
            if aa not in amino_acids:
                amino_acids.append(aa)

    if not amino_acids:
        return

    aa_totals = {aa: sum(stats_per_pos[pos]['aa_counts'].get(aa, 0) for pos in positions)
                 for aa in amino_acids}
    amino_acids.sort(key=lambda aa: aa_totals[aa], reverse=True)

    matrix = np.zeros((len(amino_acids), len(positions)), dtype=float)
    for j, pos in enumerate(positions):
        aa_counts = stats_per_pos[pos]['aa_counts']
        for i, aa in enumerate(amino_acids):
            matrix[i, j] = aa_counts.get(aa, 0)

    fig_width = max(6, len(positions) * 0.7)
    fig = plt.figure(figsize=(fig_width, 5.2))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.2, 1.2], hspace=0.25)

    ax_heat = fig.add_subplot(gs[0])
    ax_bottom = fig.add_subplot(gs[1], sharex=ax_heat)

    ax_heat.set_yticks(range(len(amino_acids)))
    ax_heat.set_yticklabels(amino_acids, fontsize=9)
    ax_heat.set_xticks(range(len(positions)))
    ax_heat.set_xticklabels(positions, rotation=45, ha='right')
    ax_heat.set_ylabel('Amino acid')
    ax_heat.set_xlim(-0.5, len(positions) - 0.5)
    ax_heat.set_ylim(-0.5, len(amino_acids) - 0.5)
    ax_heat.invert_yaxis()

    for i in range(len(amino_acids)):
        for j in range(len(positions)):
            value = matrix[i, j]
            if value > 0:
                display_value = int(value) if value < 8000 else f">{int(value//1000)}k"
                facecolor = 'white'
            else:
                display_value = ''
                facecolor = '#f6f6f6'
            rect = plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                 facecolor=facecolor, edgecolor='#b0b0b0', linewidth=0.6)
            ax_heat.add_patch(rect)
            if display_value:
                ax_heat.text(j, i, display_value, ha='center', va='center', fontsize=8,
                             family='monospace')

    ax_heat.axvline(-0.5, color='#7f7f7f', linestyle=':', linewidth=0.8)
    ax_heat.axvline(len(positions) - 0.5, color='#7f7f7f', linestyle=':', linewidth=0.8)
    ax_heat.text((len(positions) - 1) / 2, -1.1,
                 f'Window {positions[0]} – {positions[-1]}',
                 ha='center', va='top', fontsize=9)

    categories = []
    for pos in positions:
        cat = stats_per_pos[pos]['category']
        if cat not in ('synonymous-only', 'nonsynonymous-only'):
            cat = 'invariant'
        categories.append(cat)

    category_colors = {
        'synonymous-only': CATEGORY_COLORS['synonymous-only'],
        'nonsynonymous-only': CATEGORY_COLORS['nonsynonymous-only'],
        'invariant': '#9e9e9e',
    }

    n_variants = [stats_per_pos[pos]['n_variants'] for pos in positions]
    colors = [category_colors[cat] for cat in categories]
    ax_bottom.axvspan(-0.5, len(positions) - 0.5, color='#f2f2f2', alpha=0.25)
    ax_bottom.bar(range(len(positions)), n_variants, color='white', edgecolor=colors, linewidth=1.5)
    ax_bottom.set_ylabel('Distinct codons', labelpad=8)
    ax_bottom.set_xlabel('Codon position', labelpad=10)
    ax_bottom.set_xticks(range(len(positions)))
    ax_bottom.set_xticklabels(positions, rotation=45, ha='right')
    ax_bottom.axvline(focus_col, color='black', linestyle='--', linewidth=1)
    legend_handles = [
        plt.Line2D([], [], marker='s', linestyle='', color='white', markeredgecolor=category_colors[key],
                   label=key.replace('-', ' '))
        for key in ['invariant', 'synonymous-only', 'nonsynonymous-only']
    ]
    ax_bottom.legend(handles=legend_handles, frameon=False, loc='upper right', fontsize=8)

    fig.tight_layout()
    out_path = output_dir / f'{gene}_codon_{focus_idx}_window.pdf'
    fig.savefig(out_path)
    plt.close(fig)


def ensure_outputs(output_dir: Path) -> None:
    (output_dir / 'tables').mkdir(parents=True, exist_ok=True)
    (output_dir / 'plots').mkdir(parents=True, exist_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Visualise variable codon positions in operon alignments.'
    )
    parser.add_argument(
        '--alignment-dir',
        type=Path,
        default=Path('../05_operon_assembly_extraction/output/msa/dna_alignments'),
        help='Directory containing gene alignments (FASTA)'
    )
    parser.add_argument(
        '--genes',
        nargs='*',
        default=['ptsA'],
        help='Gene names to process (defaults to ptsA)'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('output'),
        help='Output directory for tables and plots'
    )
    parser.add_argument(
        '--report',
        action='store_true',
        help='Write a text report summarising highlighted codons'
    )
    parser.add_argument(
        '--focus-codons',
        type=int,
        nargs='*',
        default=None,
        help='Optional codon indices to generate per-site allele plots'
    )
    parser.add_argument(
        '--max-examples',
        type=int,
        default=3,
        help='Maximum sequences shown per codon in raw alignment examples'
    )
    parser.add_argument(
        '--window-size',
        type=int,
        default=5,
        help='Half-window size (in codons) for neighbourhood heatmaps around focus sites'
    )

    args = parser.parse_args()

    ensure_outputs(args.output_dir)
    report_lines: List[str] = []

    for gene in args.genes:
        alignment_path = args.alignment_dir / f'{gene}_aligned.fasta'
        if not alignment_path.exists():
            print(f'[warning] Alignment not found for {gene}: {alignment_path}')
            continue

        df, records, alignment = summarise_gene(alignment_path)
        if df.empty:
            print(f'[info] No variable codons detected for {gene}')
            continue

        table_path = args.output_dir / 'tables' / f'{gene}_variable_codons.csv'
        df.to_csv(table_path, index=False)
        print(f'[done] Wrote codon table for {gene} -> {table_path}')

        total_codons = alignment.get_alignment_length() // 3
        highlight_ranges = []
        if args.focus_codons:
            for codon in args.focus_codons:
                try:
                    idx = int(codon)
                except (TypeError, ValueError):
                    print(f'[warning] Skipping non-integer focus codon {codon!r} for {gene}')
                    continue
                start = max(1, idx - args.window_size)
                end = min(total_codons, idx + args.window_size)
                highlight_ranges.append((start - 0.5, end + 0.5))

        plot_gene(df, gene, args.output_dir / 'plots', highlight_ranges if highlight_ranges else None)
        plot_amino_acid_diversity(records, gene, args.output_dir / 'plots')
        print(f'[done] Wrote codon plot for {gene}')

        if args.focus_codons:
            plot_focus_codons(
                records,
                alignment,
                args.focus_codons,
                gene,
                args.output_dir / 'plots' / 'focus_sites',
                max_examples=args.max_examples,
                window_size=args.window_size
            )
            codon_dir = args.output_dir / 'plots' / 'focus_sites' / 'windows_codons'
            for codon in args.focus_codons:
                if alignment.get_alignment_length() // 3 >= codon:
                    plot_codon_window_codons(
                        alignment,
                        codon,
                        gene,
                        codon_dir,
                        window_size=args.window_size
                    )
            combined_dir = args.output_dir / 'plots' / 'focus_sites' / 'windows_combined'
            for codon in args.focus_codons:
                if alignment.get_alignment_length() // 3 >= codon:
                    plot_codon_window_combined(
                        alignment,
                        codon,
                        gene,
                        combined_dir,
                        window_size=args.window_size
                    )

        if args.report:
            n_syn = int((df['category'] == 'synonymous-only').sum())
            n_nonsyn = int((df['category'] == 'nonsynonymous-only').sum())
            top_rows = df.nlargest(5, 'n_variants')
            report_lines.append(
                f"Gene {gene}\n"
                f"  Variable codons: {len(df)}\n"
                f"  Synonymous-only: {n_syn}\n"
                f"  Non-synonymous-only: {n_nonsyn}\n"
                f"  Top sites: " + "; ".join(
                    f"{row['codon_index']} ({row['category']}, {row['aa_states']})"
                    for _, row in top_rows.iterrows()
                )
            )

    if args.report and report_lines:
        report_path = args.output_dir / 'codon_variation_highlights.txt'
        report_path.write_text("\n\n".join(report_lines))
        print(f'[done] Wrote report -> {report_path}')


if __name__ == '__main__':
    main()
