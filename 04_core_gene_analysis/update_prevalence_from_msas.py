#!/usr/bin/env python3
"""
Update gene prevalence stats using filtered MSA sequence counts.
Regenerates threshold curve and summary based on 8,596 genomes.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import os

TOTAL_GENOMES = 8596
OUTPUT_DIR = 'output'

def update_prevalence_stats():
    """Update gene_prevalence_stats.csv using MSA sequence counts."""
    print("=" * 60)
    print("Updating Gene Prevalence Stats")
    print("=" * 60)

    # Load filtered MSA metrics (has correct n_sequences)
    msa_metrics = pd.read_csv(f'{OUTPUT_DIR}/core_gene_conservation_metrics.csv')
    print(f"Loaded {len(msa_metrics)} genes from MSA metrics")

    # Load existing prevalence stats (has gene descriptions)
    prevalence_file = f'{OUTPUT_DIR}/gene_prevalence_stats.csv'
    old_prevalence = pd.read_csv(prevalence_file)
    print(f"Loaded {len(old_prevalence)} genes from existing prevalence stats")

    # Create gene name mapping (remove _aligned suffix)
    msa_metrics['gene_base'] = msa_metrics['gene'].str.replace('_aligned', '')

    # Merge to get updated counts
    merged = old_prevalence.merge(
        msa_metrics[['gene_base', 'n_sequences']],
        left_on='gene',
        right_on='gene_base',
        how='left'
    )

    # Update counts where we have MSA data
    merged['count'] = merged['n_sequences'].fillna(merged['count'])

    # For genes not in MSAs (non-core genes), estimate by subtracting 7
    # (conservative: assumes excluded genomes had all genes)
    mask = merged['n_sequences'].isna()
    merged.loc[mask, 'count'] = (merged.loc[mask, 'count'] - 7).clip(lower=0)

    # Recalculate prevalence
    merged['count'] = merged['count'].astype(int)
    merged['prevalence'] = merged['count'] / TOTAL_GENOMES

    # Keep only original columns
    result = merged[['gene', 'count', 'prevalence', 'product']].copy()
    result = result.sort_values('count', ascending=False)

    # Save updated prevalence
    result.to_csv(prevalence_file, index=False)
    print(f"\nUpdated {prevalence_file}")
    print(f"  Max count: {result['count'].max()} (was 8603)")
    print(f"  Genes with 100% prevalence: {(result['prevalence'] == 1.0).sum()}")

    return result


def update_threshold_summary(prevalence_df):
    """Regenerate threshold summary CSV."""
    print("\n" + "=" * 60)
    print("Updating Threshold Summary")
    print("=" * 60)

    thresholds = list(range(0, 101))
    results = []

    for thresh in thresholds:
        min_count = int(TOTAL_GENOMES * thresh / 100)
        gene_count = (prevalence_df['count'] >= min_count).sum()
        results.append({
            'threshold_pct': thresh,
            'gene_count': gene_count
        })

    df = pd.DataFrame(results)

    # Save
    summary_file = f'{OUTPUT_DIR}/core_gene_threshold_summary.csv'
    df.to_csv(summary_file, index=False)

    # Key stats
    core_95 = df[df['threshold_pct'] == 95]['gene_count'].values[0]
    print(f"Saved {summary_file}")
    print(f"  At 95% threshold: {core_95} core genes")

    return df


def update_core_genes_list(prevalence_df):
    """Update list of core genes at 95% threshold."""
    print("\n" + "=" * 60)
    print("Updating Core Genes List")
    print("=" * 60)

    min_count = int(TOTAL_GENOMES * 0.95)  # 8166
    core_genes = prevalence_df[prevalence_df['count'] >= min_count]['gene'].tolist()

    list_file = f'{OUTPUT_DIR}/core_genes_95pct.txt'
    with open(list_file, 'w') as f:
        for gene in sorted(core_genes):
            f.write(f"{gene}\n")

    print(f"Saved {list_file}")
    print(f"  Core genes at 95%: {len(core_genes)}")

    return core_genes


def regenerate_threshold_plot(threshold_df):
    """Regenerate the threshold curve plot."""
    print("\n" + "=" * 60)
    print("Regenerating Threshold Curve Plot")
    print("=" * 60)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Full range plot
    ax1.plot(threshold_df['threshold_pct'], threshold_df['gene_count'],
             'b-', linewidth=2)
    ax1.axvline(x=95, color='red', linestyle='--', alpha=0.7, label='95% threshold')
    ax1.axhline(y=threshold_df[threshold_df['threshold_pct']==95]['gene_count'].values[0],
                color='red', linestyle='--', alpha=0.7)
    ax1.set_xlabel('Prevalence Threshold (%)', fontsize=12)
    ax1.set_ylabel('Number of Core Genes', fontsize=12)
    ax1.set_title(f'Core Gene Count vs Prevalence Threshold\n(n={TOTAL_GENOMES:,} genomes)', fontsize=12)
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Zoomed plot (90-100%)
    zoom_df = threshold_df[threshold_df['threshold_pct'] >= 90]
    ax2.plot(zoom_df['threshold_pct'], zoom_df['gene_count'],
             'b-', linewidth=2, marker='o', markersize=6)
    ax2.axvline(x=95, color='red', linestyle='--', alpha=0.7, label='95% threshold')

    core_95 = threshold_df[threshold_df['threshold_pct']==95]['gene_count'].values[0]
    ax2.axhline(y=core_95, color='red', linestyle='--', alpha=0.7)
    ax2.annotate(f'{core_95:,} genes', xy=(95, core_95),
                 xytext=(96, core_95+20), fontsize=10,
                 arrowprops=dict(arrowstyle='->', color='red'))

    ax2.set_xlabel('Prevalence Threshold (%)', fontsize=12)
    ax2.set_ylabel('Number of Core Genes', fontsize=12)
    ax2.set_title('High Prevalence Region (90-100%)', fontsize=12)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(89, 101)

    plt.tight_layout()

    # Save
    plot_file = f'{OUTPUT_DIR}/core_gene_threshold_curve.png'
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    plt.savefig(plot_file.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()

    print(f"Saved {plot_file}")
    print(f"Saved {plot_file.replace('.png', '.pdf')}")


def main():
    os.chdir('/work')

    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Target genome count: {TOTAL_GENOMES:,}")

    # Step 1: Update prevalence stats
    prevalence_df = update_prevalence_stats()

    # Step 2: Update threshold summary
    threshold_df = update_threshold_summary(prevalence_df)

    # Step 3: Update core genes list
    core_genes = update_core_genes_list(prevalence_df)

    # Step 4: Regenerate plot
    regenerate_threshold_plot(threshold_df)

    print("\n" + "=" * 60)
    print("DONE!")
    print("=" * 60)


if __name__ == '__main__':
    main()
