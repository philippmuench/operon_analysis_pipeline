#!/usr/bin/env python3
"""Plot mismatch/deletion windows around known variation sites."""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List, Set

import matplotlib.pyplot as plt
import pandas as pd


def load_summary(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = {
        "position",
        "mismatch_rate_legacy",
        "mismatch_rate_updated",
        "deletion_rate_legacy",
        "deletion_rate_updated",
        "coverage_legacy",
        "coverage_updated",
        "known_variation",
        "feature_gene",
        "feature_codon",
    }
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Summary {path} missing columns: {sorted(missing)}")
    return df


def collect_variation_positions(df: pd.DataFrame, specific: Iterable[int] | None = None) -> List[int]:
    if specific:
        return sorted(set(int(pos) for pos in specific))
    mask = df["known_variation"].fillna("").astype(str).str.len() > 0
    return sorted(int(pos) for pos in df.loc[mask, "position"].unique())


def plot_window(
    df: pd.DataFrame,
    center: int,
    window: int,
    out_dir: Path,
    label: str,
    legacy_label: str,
    updated_label: str,
) -> None:
    start = center - window
    end = center + window
    window_df = df[(df["position"] >= start) & (df["position"] <= end)].copy()
    if window_df.empty:
        return
    window_df.sort_values("position", inplace=True)
    window_df["mismatch_rate_legacy_pct"] = window_df["mismatch_rate_legacy"] * 100
    window_df["mismatch_rate_updated_pct"] = window_df["mismatch_rate_updated"] * 100
    window_df["deletion_rate_legacy_pct"] = window_df["deletion_rate_legacy"] * 100
    window_df["deletion_rate_updated_pct"] = window_df["deletion_rate_updated"] * 100

    center_row = window_df.loc[window_df["position"] == center].iloc[0]
    gene = center_row.get("feature_gene", "")
    codon = center_row.get("feature_codon", "")
    var_label = center_row.get("known_variation", "")

    fig, axes = plt.subplots(2, 1, figsize=(8, 5), sharex=True)

    axes[0].plot(
        window_df["position"],
        window_df["mismatch_rate_legacy_pct"],
        label=f"Mismatch {legacy_label}",
        color="#4C72B0",
    )
    axes[0].plot(
        window_df["position"],
        window_df["mismatch_rate_updated_pct"],
        label=f"Mismatch {updated_label}",
        color="#DD8452",
    )
    axes[0].axvline(center, color="black", linestyle="--", linewidth=1)
    axes[0].set_ylabel("Mismatch rate (%)")
    axes[0].legend(loc="upper right")

    axes[1].plot(
        window_df["position"],
        window_df["deletion_rate_legacy_pct"],
        label=f"Deletion {legacy_label}",
        color="#55A868",
    )
    axes[1].plot(
        window_df["position"],
        window_df["deletion_rate_updated_pct"],
        label=f"Deletion {updated_label}",
        color="#C44E52",
    )
    axes[1].axvline(center, color="black", linestyle="--", linewidth=1)
    axes[1].set_ylabel("Deletion rate (%)")
    axes[1].set_xlabel("Operon position (nt)")
    axes[1].legend(loc="upper right")

    title_parts = [f"Position {center}"]
    if gene and isinstance(gene, str) and gene.strip():
        title_parts.append(gene)
    if not pd.isna(codon):
        title_parts.append(f"codon {int(codon)}")
    if var_label:
        title_parts.append(f"variant {var_label}")

    fig.suptitle(" | ".join(title_parts))
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    filename = f"variation_{center}_{gene or 'unknown'}.png"
    filename = filename.replace("/", "_")
    fig.savefig(out_dir / filename, dpi=200)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--combined",
        type=Path,
        default=Path("output_full_run/position_stats/combined_position_summary.tsv"),
        help="Annotated combined position summary TSV",
    )
    parser.add_argument(
        "--positions",
        type=int,
        nargs="*",
        help="Specific positions to plot (default: use rows with known_variation)",
    )
    parser.add_argument(
        "--window",
        type=int,
        default=10,
        help="Number of nucleotides to include on each side (default: 10)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output_full_run/position_stats/variation_windows"),
        help="Directory to write plots",
    )
    parser.add_argument(
        "--legacy-label",
        default="Legacy",
        help="Label for legacy dataset",
    )
    parser.add_argument(
        "--updated-label",
        default="Updated",
        help="Label for updated dataset",
    )
    args = parser.parse_args()

    df = load_summary(args.combined)
    positions = collect_variation_positions(df, args.positions)
    if not positions:
        raise SystemExit("No variation positions found")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    for pos in positions:
        plot_window(
            df,
            pos,
            args.window,
            args.output_dir,
            label=str(pos),
            legacy_label=args.legacy_label,
            updated_label=args.updated_label,
        )


if __name__ == "__main__":
    main()
