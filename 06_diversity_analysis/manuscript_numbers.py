#!/usr/bin/env python3
"""
Generate manuscript statistics for comparative diversity analysis.
"""

import os
import pandas as pd
from datetime import datetime


def generate_manuscript_stats():
    """Generate statistics for manuscript."""

    print("Comparative Diversity Analysis Statistics for Manuscript")
    print("=" * 60)
    print(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Load core gene metrics
    core_csv = "../04_core_gene_analysis/output/core_gene_conservation_metrics.csv"
    operon_csv = "output/operon_conservation_summary.csv"

    print(f"\n1. Input Data:")

    if os.path.exists(core_csv):
        core_df = pd.read_csv(core_csv)
        print(f"   Core genes analyzed: {len(core_df):,}")
        print(f"   Mean core gene conservation: {core_df['mean_conservation'].mean():.3f}")
        print(f"   Core gene conservation range: {core_df['mean_conservation'].min():.3f} - {core_df['mean_conservation'].max():.3f}")
    else:
        print(f"   WARNING: Core gene metrics not found at {core_csv}")

    print(f"\n2. Operon Gene Conservation Rankings:")

    if os.path.exists(operon_csv):
        operon_df = pd.read_csv(operon_csv)
        print(f"   Operon genes analyzed: {len(operon_df)}")
        print(f"   Total core genes for ranking: {operon_df['Total Core Genes'].iloc[0]:,}")

        print(f"\n   Individual gene rankings:")
        for _, row in operon_df.iterrows():
            print(f"     {row['Gene']}: {row['Conservation (%)']:.2f}% conservation, Rank {row['Rank']}/{row['Total Core Genes']} ({row['Percentile']:.1f}th percentile)")

        # Summary statistics
        print(f"\n3. Summary Statistics:")
        print(f"   Mean operon conservation: {operon_df['Conservation (%)'].mean():.2f}%")
        print(f"   Median operon conservation: {operon_df['Conservation (%)'].median():.2f}%")
        print(f"   Conservation range: {operon_df['Conservation (%)'].min():.2f}% - {operon_df['Conservation (%)'].max():.2f}%")
        print(f"   Mean percentile rank: {operon_df['Percentile'].mean():.1f}")

        # Most and least conserved
        most_conserved = operon_df.loc[operon_df['Conservation (%)'].idxmax()]
        least_conserved = operon_df.loc[operon_df['Conservation (%)'].idxmin()]
        print(f"\n   Most conserved operon gene: {most_conserved['Gene']} ({most_conserved['Conservation (%)']:.2f}%, rank {most_conserved['Rank']})")
        print(f"   Least conserved operon gene: {least_conserved['Gene']} ({least_conserved['Conservation (%)']:.2f}%, rank {least_conserved['Rank']})")

        # Compare to core gene average
        if os.path.exists(core_csv):
            core_mean = core_df['mean_conservation'].mean() * 100  # Convert to percentage
            operon_mean = operon_df['Conservation (%)'].mean()
            diff = operon_mean - core_mean
            print(f"\n4. Comparison to Core Genes:")
            print(f"   Core gene mean conservation: {core_mean:.2f}%")
            print(f"   Operon gene mean conservation: {operon_mean:.2f}%")
            print(f"   Difference: {diff:+.2f}%")
            if diff > 0:
                print(f"   Operon genes are {abs(diff):.2f}% MORE conserved than average core genes")
            else:
                print(f"   Operon genes are {abs(diff):.2f}% LESS conserved than average core genes")
    else:
        print(f"   WARNING: Operon summary not found at {operon_csv}")

    print(f"\n5. Output Files:")
    print(f"   conservation_ranking_bars.png/pdf - Bar plot of conservation rankings")
    print(f"   conservation_distribution.png/pdf - Distribution histogram")
    print(f"   operon_conservation_summary.csv - Summary table")

    print("\n" + "=" * 60)
    print("MANUSCRIPT NUMBERS SUMMARY:")
    print("=" * 60)

    if os.path.exists(operon_csv):
        operon_df = pd.read_csv(operon_csv)
        total_core = operon_df['Total Core Genes'].iloc[0]
        print(f"Total core genes analyzed: {total_core:,}")
        print(f"Operon genes analyzed: {len(operon_df)}")
        print(f"Mean operon conservation: {operon_df['Conservation (%)'].mean():.2f}%")
        print(f"Operon genes are more conserved than {operon_df['Percentile'].mean():.1f}% of core genes on average")

        if os.path.exists(core_csv):
            core_mean = core_df['mean_conservation'].mean() * 100
            operon_mean = operon_df['Conservation (%)'].mean()
            diff = operon_mean - core_mean
            print(f"Operon vs core gene conservation difference: {diff:+.2f}%")


if __name__ == "__main__":
    import sys

    # Change to script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)

    # Redirect output to file if specified
    if len(sys.argv) > 1:
        output_file = sys.argv[1]
        with open(output_file, 'w') as f:
            import contextlib
            with contextlib.redirect_stdout(f):
                generate_manuscript_stats()
        print(f"Results saved to: {output_file}")
    else:
        generate_manuscript_stats()
