#!/usr/bin/env python3
"""Visualise full-operon BLAST mapping results."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def identity_scatter(df: pd.DataFrame, out_path: Path) -> None:
    plt.figure(figsize=(6, 6))
    plt.scatter(df["legacy_pident"], df["updated_pident"], alpha=0.7, s=10)
    min_val = min(df["legacy_pident"].min(), df["updated_pident"].min())
    max_val = max(df["legacy_pident"].max(), df["updated_pident"].max())
    plt.plot([min_val, max_val], [min_val, max_val], linestyle="--", color="grey", linewidth=1)
    plt.xlabel("Legacy identity (%)")
    plt.ylabel("Updated identity (%)")
    plt.title("Full operon identity: updated vs legacy")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mapping_tsv", type=Path, help="Path to full_operon_mapping.tsv")
    parser.add_argument("--output-dir", type=Path, default=Path("plots"))
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.mapping_tsv, sep="\t")
    df["legacy_found"] = df["legacy_found"].fillna(False).astype(bool)
    df["updated_found"] = df["updated_found"].fillna(False).astype(bool)

    both_hits = df[df["legacy_found"] & df["updated_found"]]
    identity_scatter(both_hits, args.output_dir / "identity_scatter.png")


if __name__ == "__main__":
    main()
