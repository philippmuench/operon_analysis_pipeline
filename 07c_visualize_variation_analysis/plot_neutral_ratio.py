#!/usr/bin/env python3
"""Visualise why the neutral nonsyn/syn ratio is ~2.5."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec


GENETIC_CODE = {
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
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}


def codon_opportunities() -> tuple[int, int, list[float]]:
    bases = 'ACGT'
    syn_total = 0
    nonsyn_total = 0
    ratios: list[float] = []

    for codon, aa in GENETIC_CODE.items():
        if aa == '*':
            continue
        syn = 0
        nonsyn = 0
        for pos in range(3):
            for base in bases:
                if base == codon[pos]:
                    continue
                mutated = codon[:pos] + base + codon[pos + 1:]
                aa_mut = GENETIC_CODE.get(mutated)
                if aa_mut is None or aa_mut == '*':
                    continue
                if aa_mut == aa:
                    syn += 1
                else:
                    nonsyn += 1
        syn_total += syn
        nonsyn_total += nonsyn
        if syn > 0:
            ratios.append(nonsyn / syn)

    return syn_total, nonsyn_total, ratios


def neighbour_summary(codon: str, bases: str = 'ACGT'):
    aa = GENETIC_CODE[codon]
    entries = []
    syn = nonsyn = 0
    for pos in range(3):
        alt_bases = [b for b in bases if b != codon[pos]]
        for row, base in enumerate(alt_bases):
            mutated = codon[:pos] + base + codon[pos + 1:]
            aa_mut = GENETIC_CODE.get(mutated)
            if aa_mut is None or aa_mut == '*':
                kind = 'stop'
            elif aa_mut == aa:
                kind = 'syn'
                syn += 1
            else:
                kind = 'nonsyn'
                nonsyn += 1
            entries.append((pos, row, mutated, aa_mut, kind))
    return syn, nonsyn, entries


def plot_codon_grid(ax, codon: str) -> None:
    syn, nonsyn, entries = neighbour_summary(codon)
    aa = GENETIC_CODE[codon]
    colors = {'syn': '#4daf4a', 'nonsyn': '#e41a1c', 'stop': '#c7c7c7'}

    ax.set_xlim(-0.5, 2.5)
    ax.set_ylim(-0.5, 2.5)
    ax.invert_yaxis()
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(['Position 1', 'Position 2', 'Position 3'], fontsize=8)
    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels(['Alt 1', 'Alt 2', 'Alt 3'], fontsize=8)
    ax.tick_params(length=0)

    for pos, row, mutated, aa_mut, kind in entries:
        rect = plt.Rectangle((pos - 0.5, row - 0.5), 1, 1,
                             facecolor=colors[kind], edgecolor='black', linewidth=0.8,
                             alpha=0.75 if kind != 'stop' else 0.25)
        ax.add_patch(rect)
        label = mutated if kind != 'stop' else f"{mutated}\n(stop)"
        aa_label = aa_mut if aa_mut and aa_mut != '*' else ''
        ax.text(pos, row, f"{label}\n({aa_label})" if aa_label else label,
                ha='center', va='center', fontsize=8, family='monospace')

    ax.set_title(f"{codon} → {aa}\nSyn={syn}, Nonsyn={nonsyn}", fontsize=9, fontweight='bold')


def main() -> None:
    syn_total, nonsyn_total, ratios = codon_opportunities()
    neutral_ratio = nonsyn_total / syn_total

    out_dir = Path('output/plots')
    out_dir.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(14, 4))
    outer = GridSpec(1, 3, width_ratios=[1.1, 1.8, 1.3], wspace=0.35)

    ax_bar = fig.add_subplot(outer[0])
    bars = ax_bar.bar(['Synonymous', 'Non-synonymous'], [syn_total, nonsyn_total],
                      color=['#4daf4a', '#e41a1c'])
    ax_bar.set_ylabel('Single-nucleotide opportunities')
    ax_bar.set_title('All sense codons')
    for bar, value in zip(bars, [syn_total, nonsyn_total]):
        ax_bar.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + max(syn_total, nonsyn_total) * 0.025,
                    f'{value}', ha='center', va='bottom', fontsize=10, fontweight='bold')
    ax_bar.text(0.5, max(syn_total, nonsyn_total) * 1.1,
                f'Neutral ratio = {neutral_ratio:.2f}',
                ha='center', fontsize=11, fontweight='bold')

    middle = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[1], wspace=0.3)
    ax_cod1 = fig.add_subplot(middle[0])
    ax_cod2 = fig.add_subplot(middle[1])
    plot_codon_grid(ax_cod1, 'GCT')  # Ala, fourfold degenerate
    plot_codon_grid(ax_cod2, 'ATG')  # Met, no synonymous neighbours

    ax_hist = fig.add_subplot(outer[2])
    ax_hist.hist(ratios, bins=15, color='#377eb8', edgecolor='black', alpha=0.7)
    ax_hist.set_xlabel('Codon nonsyn/syn opportunity ratio')
    ax_hist.set_ylabel('Number of codons')
    ax_hist.set_title('Distribution across 61 sense codons')
    ax_hist.axvline(neutral_ratio, color='black', linestyle='--', linewidth=1,
                    label=f'Mean = {neutral_ratio:.2f}')
    ax_hist.legend(frameon=False, loc='upper right')

    fig.tight_layout()

    out_path = out_dir / 'neutral_ratio_overview.pdf'
    fig.savefig(out_path, bbox_inches='tight')
    plt.close(fig)

    print('Total synonymous opportunities:', syn_total)
    print('Total non-synonymous opportunities:', nonsyn_total)
    print('Neutral nonsyn/syn ratio:', f'{neutral_ratio:.4f}')
    print('Codons with no synonymous neighbours:', 61 - len(ratios))
    print('Plot saved to:', out_path)


if __name__ == '__main__':
    main()
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
