#!/usr/bin/env python3
"""Create per-position base-call slices (match/mismatch/gap) for key variation sites."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

DEFAULT_POSITIONS = [859, 2155, 3455, 3568, 4652, 4653, 4659, 4673, 7501, 7502]


def load_tables(summary_path: Path, counts_path: Path) -> pd.DataFrame:
    summary = pd.read_csv(summary_path, sep="\t")
    counts = pd.read_csv(counts_path, sep="\t")
    merged = summary.merge(counts, on="position", how="left", suffixes=('', '_count'))
    return merged


def plot_slice(row: pd.Series, out_dir: Path) -> None:
    pos = int(row['position'])
    gene = row.get('feature_gene', '')
    codon = row.get('feature_codon')
    var_label = row.get('known_variation', '')

    match = int(row.get('match', 0) or 0)
    mismatch = int(row.get('mismatch', 0) or 0)
    gap = int(row.get('gap', 0) or 0)
    coverage = int(row.get('coverage', 0) or 0)

    values = [match, mismatch, gap]
    labels = ['Match', 'Mismatch', 'Gap']
    colors = ['#4C72B0', '#DD8452', '#55A868']

    fig, ax = plt.subplots(figsize=(4, 4))
    if coverage > 0:
        wedges, texts, autotexts = ax.pie(values, labels=labels, colors=colors, autopct='%1.2f%%', startangle=90)
        for text in (texts + autotexts):
            text.set_fontsize(8)
    else:
        ax.pie([1], labels=['no coverage'], colors=['lightgray'])

    title_parts = [f"Pos {pos}"]
    if isinstance(gene, str) and gene.strip():
        title_parts.append(gene)
    if pd.notna(codon):
        title_parts.append(f"codon {int(codon)}")
    if isinstance(var_label, str) and var_label:
        title_parts.append(f"variant {var_label}")
    title_parts.append(f"n={coverage}")

    ax.set_title(" | ".join(title_parts), fontsize=9)
    plt.tight_layout()

    filename = f"variation_slice_{pos}_{gene or 'unknown'}.png".replace('/', '_')
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / filename, dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary', type=Path, default=Path('output_full_run/position_stats/combined_position_summary.tsv'))
    parser.add_argument('--counts', type=Path, default=Path('output_full_run/position_stats/coverage_counts.tsv'))
    parser.add_argument('--positions', type=int, nargs='*', default=None)
    parser.add_argument('--output-dir', type=Path, default=Path('output_full_run/position_stats/variation_slices'))
    args = parser.parse_args()

    merged = load_tables(args.summary, args.counts)
    positions = args.positions if args.positions else DEFAULT_POSITIONS

    for pos in positions:
        row = merged.loc[merged['position'] == pos]
        if row.empty:
            continue
        plot_slice(row.iloc[0], args.output_dir)


if __name__ == '__main__':
    main()
