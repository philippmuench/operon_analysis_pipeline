#!/usr/bin/env python3
"""Visualize codon-level variation for operon alignments.

This helper script reproduces the simple codon counting logic used in the
selection analysis and emits per-codon summaries plus an annotated plot that
highlights which positions are variable and how they vary (synonymous vs
non-synonymous). The goal is to provide ready-to-plot material for figures.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt
import pandas as pd
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import numpy as np

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
    'mixed': '#6A1B9A',
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

        if syn and nonsyn:
            category = 'mixed'
        elif syn:
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


def plot_gene(df: pd.DataFrame, gene: str, output_dir: Path) -> None:
    if df.empty:
        return

    fig, ax = plt.subplots(figsize=(14, 3.0))

    y = [1] * len(df)
    colors = [CATEGORY_COLORS.get(cat, '#555555') for cat in df['category']]
    sizes = [60 + 40 * max(0, n - 2) for n in df['n_variants']]  # Highlight diverse codons

    ax.scatter(df['codon_index'], y, c=colors, s=sizes, linewidths=0.6,
               edgecolors='black', alpha=0.9)

    # Label sites where more than two unique codons occur
    for _, row in df[df['n_variants'] > 2].iterrows():
        ax.text(row['codon_index'], 1, str(int(row['n_variants'])),
                fontsize=8, fontweight='bold', color='white', ha='center', va='center')

    ax.set_title(f'{gene} – Variable Codon Positions', fontsize=14, fontweight='bold')
    ax.set_xlabel('Codon position')
    ax.set_yticks([])
    ax.set_xlim(df['codon_index'].min() - 1, df['codon_index'].max() + 1)
    ax.grid(axis='x', alpha=0.2)

    # Annotate up to five largest non-synonymous sites for storytelling
    nonsyn = df[df['nonsyn_changes'] == 1]
    for _, row in nonsyn.nlargest(5, 'n_variants').iterrows():
        ax.text(row['codon_index'], 1.08,
                f"AA: {row['aa_states']}",
                fontsize=8, ha='center', va='bottom', rotation=0)

    # Legend keyed by category
    handles = [
        plt.Line2D([], [], marker='o', linestyle='', color=color,
                   markeredgecolor='black', label=label.replace('-', ' '))
        for label, color in CATEGORY_COLORS.items()
    ]
    ax.legend(handles=handles, loc='upper right', frameon=False)

    fig.tight_layout()
    out_path = output_dir / f'{gene}_variable_codons.pdf'
    fig.savefig(out_path)
    plt.close(fig)


def plot_focus_codons(records: List[CodonSummary], alignment: MultipleSeqAlignment,
                      focus_codons: Iterable[int], gene: str, output_dir: Path,
                      max_examples: int = 3, window_size: int = 5) -> None:
    if not focus_codons:
        return

    record_lookup = {rec.codon_index: rec for rec in records}
    output_dir.mkdir(parents=True, exist_ok=True)
    example_dir = output_dir / 'examples'
    example_dir.mkdir(parents=True, exist_ok=True)
    window_dir = output_dir / 'windows'
    window_dir.mkdir(parents=True, exist_ok=True)
    palette = plt.get_cmap('tab20')

    for codon_idx in focus_codons:
        rec = record_lookup.get(int(codon_idx))
        if rec is None:
            print(f'[warning] Focus codon {codon_idx} not found for {gene}')
            continue

        items = sorted(rec.codon_counts.items(), key=lambda kv: (-kv[1], kv[0]))
        aas = [GENETIC_CODE[codon] for codon, _ in items]
        unique_aas = list(dict.fromkeys(aas))
        aa_colors = {aa: palette(i % palette.N) for i, aa in enumerate(unique_aas)}

        counts = [count for _, count in items]
        labels = [f"{codon}\n({GENETIC_CODE[codon]})" for codon, _ in items]
        colors = [aa_colors[aa] for aa in aas]

        fig, ax = plt.subplots(figsize=(6, 4))
        bars = ax.bar(range(len(items)), counts, color=colors, edgecolor='black')

        for bar, count in zip(bars, counts):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + max(counts) * 0.02,
                    str(count), ha='center', va='bottom', fontsize=9,
                    fontweight='bold')

        ax.set_xticks(range(len(items)))
        ax.set_xticklabels(labels)
        ax.set_ylabel('Sequence count')
        ax.set_title(
            f'{gene} codon {codon_idx} – {rec.category} (variants: {rec.n_variants})',
            fontweight='bold'
        )

        legend_handles = [
            plt.Line2D([], [], marker='s', linestyle='', color=color,
                       label=f'{aa} (aa)')
            for aa, color in aa_colors.items()
        ]
        ax.legend(handles=legend_handles, frameon=False, loc='upper right')
        ax.margins(y=0.1)
        fig.tight_layout()

        out_path = output_dir / f'{gene}_codon_{codon_idx}_profile.pdf'
        fig.savefig(out_path)
        plt.close(fig)

        plot_alignment_examples(alignment, codon_idx, items, aa_colors, gene,
                                 example_dir, max_examples)
        plot_codon_window_heatmap(alignment, codon_idx, gene, window_dir,
                                   window_size)


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


def compute_codon_stats(alignment: MultipleSeqAlignment, codon_idx: int) -> Dict[str, object]:
    start = (codon_idx - 1) * 3
    end = start + 3

    if start < 0 or end > alignment.get_alignment_length():
        return {
            'codon_counts': Counter(),
            'aa_counts': Counter(),
            'n_variants': 0,
            'category': 'invariant',
        }

    codon_counts: Counter[str] = Counter()
    for record in alignment:
        codon = str(record.seq[start:end]).upper()
        if len(codon) != 3 or '-' in codon or 'N' in codon or codon not in GENETIC_CODE:
            continue
        codon_counts[codon] += 1

    aa_counts: Counter[str] = Counter()
    for codon, count in codon_counts.items():
        aa_counts[GENETIC_CODE[codon]] += count

    n_variants = len(codon_counts)
    category = 'invariant'
    if n_variants >= 2:
        n_amino_acids = len(aa_counts)
        if n_amino_acids == 1:
            category = 'synonymous-only'
        else:
            category = 'nonsynonymous-only'

    return {
        'codon_counts': codon_counts,
        'aa_counts': aa_counts,
        'n_variants': max(n_variants, 1 if codon_counts else 0),
        'category': category,
    }


def plot_codon_window_heatmap(alignment: MultipleSeqAlignment, focus_idx: int,
                              gene: str, output_dir: Path, window_size: int) -> None:
    total_codons = alignment.get_alignment_length() // 3
    window_min = max(1, focus_idx - window_size)
    window_max = min(total_codons, focus_idx + window_size)
    positions = list(range(window_min, window_max + 1))

    stats_per_pos = {pos: compute_codon_stats(alignment, pos) for pos in positions}

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
    fig = plt.figure(figsize=(fig_width, 5.6))
    gs = fig.add_gridspec(3, 1, height_ratios=[0.7, 3.2, 1.2], hspace=0.25)

    ax_top = fig.add_subplot(gs[0])
    ax_heat = fig.add_subplot(gs[1], sharex=ax_top)
    ax_bottom = fig.add_subplot(gs[2], sharex=ax_top)

    ax_top.set_xlim(-0.5, len(positions) - 0.5)
    ax_top.set_ylim(0, 1)
    ax_top.axis('off')
    for j, pos in enumerate(positions):
        codon_counts = stats_per_pos[pos]['codon_counts']
        if not codon_counts:
            text = '—'
        else:
            top_items = sorted(codon_counts.items(), key=lambda kv: (-kv[1], kv[0]))[:2]
            formatted = []
            for cod, cnt in top_items:
                formatted.append(f"{cod} (>{cnt//1000}k)" if cnt >= 8000 else f"{cod} ({int(cnt)})")
            text = '\n'.join(formatted)
        ax_top.text(j, 0.5, text, ha='center', va='center', fontsize=7, family='monospace')

    ax_heat.set_yticks(range(len(amino_acids)))
    ax_heat.set_yticklabels(amino_acids, fontsize=9)
    ax_heat.set_xticks(range(len(positions)))
    ax_heat.set_xticklabels(positions, rotation=45, ha='right')
    ax_heat.set_ylabel('Amino acid')
    ax_heat.set_title(f'{gene} codon window around position {focus_idx}', fontsize=12)
    ax_heat.set_xlim(-0.5, len(positions) - 0.5)
    ax_heat.set_ylim(-0.5, len(amino_acids) - 0.5)
    ax_heat.invert_yaxis()

    for i in range(len(amino_acids)):
        for j in range(len(positions)):
            value = matrix[i, j]
            if value <= 0:
                continue
            display_value = int(value) if value < 8000 else f">{int(value//1000)}k"
            rect = plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                 facecolor='white', edgecolor='black', linewidth=0.8)
            ax_heat.add_patch(rect)
            ax_heat.text(j, i, display_value, ha='center', va='center', fontsize=8, family='monospace')

    focus_column = positions.index(focus_idx)
    ax_heat.add_patch(plt.Rectangle((focus_column - 0.5, -0.5), 1, len(amino_acids),
                                    fill=False, edgecolor='black', linewidth=1.5, linestyle='--'))

    categories = []
    for pos in positions:
        cat = stats_per_pos[pos]['category']
        if cat not in ('synonymous-only', 'nonsynonymous-only', 'mixed'):
            cat = 'invariant'
        categories.append(cat)

    category_colors = {
        'synonymous-only': CATEGORY_COLORS['synonymous-only'],
        'nonsynonymous-only': CATEGORY_COLORS['nonsynonymous-only'],
        'mixed': CATEGORY_COLORS['mixed'],
        'invariant': '#9e9e9e',
    }

    n_variants = [stats_per_pos[pos]['n_variants'] for pos in positions]
    colors = [category_colors[cat] for cat in categories]
    ax_bottom.bar(range(len(positions)), n_variants, color='white', edgecolor=colors, linewidth=1.5)
    ax_bottom.set_ylabel('Distinct codons', labelpad=8)
    ax_bottom.set_xlabel('Codon position', labelpad=10)
    ax_bottom.set_xticks(range(len(positions)))
    ax_bottom.set_xticklabels(positions, rotation=45, ha='right')
    ax_bottom.axvline(focus_column, color='black', linestyle='--', linewidth=1)

    legend_handles = [
        plt.Line2D([], [], marker='s', linestyle='', color='white', markeredgecolor=category_colors[key],
                   label=key.replace('-', ' '))
        for key in ['invariant', 'synonymous-only', 'nonsynonymous-only', 'mixed']
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

        plot_gene(df, gene, args.output_dir / 'plots')
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

        if args.report:
            n_syn = int((df['category'] == 'synonymous-only').sum())
            n_nonsyn = int((df['category'] == 'nonsynonymous-only').sum())
            n_mixed = int((df['category'] == 'mixed').sum())
            top_rows = df.nlargest(5, 'n_variants')
            report_lines.append(
                f"Gene {gene}\n"
                f"  Variable codons: {len(df)}\n"
                f"  Synonymous-only: {n_syn}\n"
                f"  Non-synonymous-only: {n_nonsyn}\n"
                f"  Mixed: {n_mixed}\n"
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
