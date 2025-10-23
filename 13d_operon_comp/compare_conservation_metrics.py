#!/usr/bin/env python3
"""Compare conservation metrics between legacy and SNP-aware operon analyses."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import pandas as pd

SUMMARY_METRICS = [
    "num_sequences",
    "alignment_length",
    "conservation_score",
    "gap_percentage",
    "pairwise_identity",
]

DETAILED_METRICS = [
    "mean_conservation",
    "median_conservation",
    "highly_conserved_pct",
    "moderately_conserved_pct",
    "variable_pct",
    "positions_with_gaps_pct",
]


def read_csv(path: Path, gene_col: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    if gene_col not in df.columns:
        raise ValueError(f"Expected column '{gene_col}' in {path}")
    df[gene_col] = df[gene_col].str.lower()
    return df


def prepare_dataframe(df: pd.DataFrame, suffix: str, gene_col: str = "gene") -> pd.DataFrame:
    renamed = df.copy()
    for col in renamed.columns:
        if col == gene_col:
            continue
        renamed.rename(columns={col: f"{col}_{suffix}"}, inplace=True)
    return renamed


def merge_datasets(
    legacy: pd.DataFrame,
    updated: pd.DataFrame,
    metrics: Iterable[str],
    suffix_legacy: str,
    suffix_updated: str,
) -> pd.DataFrame:
    merged = legacy.merge(updated, on="gene", how="outer")
    for metric in metrics:
        col_a = f"{metric}_{suffix_legacy}"
        col_b = f"{metric}_{suffix_updated}"
        if col_a in merged.columns and col_b in merged.columns:
            merged[f"{metric}_difference"] = merged[col_b] - merged[col_a]
    return merged


def scatter_plot(df: pd.DataFrame, metric: str, suffix_legacy: str, suffix_updated: str, output: Path) -> None:
    col_a = f"{metric}_{suffix_legacy}"
    col_b = f"{metric}_{suffix_updated}"
    if col_a not in df or col_b not in df:
        return

    plt.figure(figsize=(6, 6))
    plt.scatter(df[col_a], df[col_b])
    max_val = max(df[col_a].max(), df[col_b].max())
    min_val = min(df[col_a].min(), df[col_b].min())
    plt.plot([min_val, max_val], [min_val, max_val], color="gray", linestyle="--", linewidth=1)

    for _, row in df.iterrows():
        plt.annotate(row["gene"], (row[col_a], row[col_b]), textcoords="offset points", xytext=(4, 4))

    plt.title(f"{metric.replace('_', ' ').title()} (updated vs legacy)")
    plt.xlabel(f"Legacy ({suffix_legacy})")
    plt.ylabel(f"Updated ({suffix_updated})")
    plt.tight_layout()
    plt.savefig(output)
    plt.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--legacy-summary",
        type=Path,
        default=Path("../05_operon_assembly_extraction/output/operon_conservation_metrics.csv"),
        help="Legacy summary metrics CSV",
    )
    parser.add_argument(
        "--legacy-detailed",
        type=Path,
        default=Path("../05_operon_assembly_extraction/output/operon_conservation_metrics_detailed.csv"),
        help="Legacy detailed metrics CSV",
    )
    parser.add_argument(
        "--updated-summary",
        type=Path,
        default=Path("../13c_new_operon_assembly_extraction/output/operon_conservation_metrics.csv"),
        help="Updated summary metrics CSV",
    )
    parser.add_argument(
        "--updated-detailed",
        type=Path,
        default=Path("../13c_new_operon_assembly_extraction/output/operon_conservation_metrics_detailed.csv"),
        help="Updated detailed metrics CSV",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output/metrics_comparison"),
        help="Directory for comparison tables and plots",
    )
    parser.add_argument("--legacy-label", default="legacy", help="Label for legacy run")
    parser.add_argument("--updated-label", default="updated", help="Label for updated run")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plots_dir = args.output_dir / "plots"
    plots_dir.mkdir(exist_ok=True)

    legacy_summary = prepare_dataframe(read_csv(args.legacy_summary, "gene"), args.legacy_label)
    updated_summary = prepare_dataframe(read_csv(args.updated_summary, "gene"), args.updated_label)

    legacy_detailed = prepare_dataframe(read_csv(args.legacy_detailed, "gene"), args.legacy_label)
    updated_detailed = prepare_dataframe(read_csv(args.updated_detailed, "gene"), args.updated_label)

    summary_comparison = merge_datasets(legacy_summary, updated_summary, SUMMARY_METRICS, args.legacy_label, args.updated_label)
    detailed_comparison = merge_datasets(legacy_detailed, updated_detailed, DETAILED_METRICS, args.legacy_label, args.updated_label)

    summary_path = args.output_dir / "summary_metrics_comparison.tsv"
    detailed_path = args.output_dir / "detailed_metrics_comparison.tsv"
    summary_comparison.to_csv(summary_path, sep="\t", index=False)
    detailed_comparison.to_csv(detailed_path, sep="\t", index=False)

    for metric in SUMMARY_METRICS:
        plot_path = plots_dir / f"summary_{metric}.png"
        scatter_plot(summary_comparison, metric, args.legacy_label, args.updated_label, plot_path)

    for metric in DETAILED_METRICS:
        plot_path = plots_dir / f"detailed_{metric}.png"
        scatter_plot(detailed_comparison, metric, args.legacy_label, args.updated_label, plot_path)

    print(f"Summary comparison written to {summary_path}")
    print(f"Detailed comparison written to {detailed_path}")
    print(f"Plots saved under {plots_dir}")


if __name__ == "__main__":
    main()
