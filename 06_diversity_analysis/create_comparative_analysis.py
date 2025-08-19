#!/usr/bin/env python3
"""
Simplified diversity analysis that reads CSV files from steps 04 and 05.
Creates comparative plots and statistical analysis without redundant calculations.
"""

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def load_core_gene_metrics(core_csv_path):
    """Load core gene metrics from step 04."""
    try:
        df = pd.read_csv(core_csv_path)
        print(f"✅ Loaded {len(df)} core genes from {core_csv_path}")
        return df
    except Exception as e:
        print(f"❌ Error loading core genes: {e}")
        return pd.DataFrame()

def load_operon_metrics(operon_csv_path, strategy_name):
    """Load operon metrics from step 05."""
    try:
        df = pd.read_csv(operon_csv_path)
        df['strategy'] = strategy_name
        print(f"✅ Loaded {len(df)} operon genes from {strategy_name}: {operon_csv_path}")
        return df
    except Exception as e:
        print(f"❌ Error loading {strategy_name} operon genes: {e}")
        return pd.DataFrame()

def create_conservation_comparison_plot(core_df, operon_dfs, output_file):
    """Create comparative conservation plot."""
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Plot 1: Conservation distribution
    ax1 = axes[0, 0]
    if not core_df.empty:
        ax1.hist(core_df['mean_conservation'], bins=50, alpha=0.7, 
                label=f'Core genes (n={len(core_df)})', color='blue')
    
    colors = ['red', 'green', 'orange', 'purple']
    for i, (strategy, df) in enumerate(operon_dfs.items()):
        if not df.empty:
            ax1.hist(df['mean_conservation'], bins=50, alpha=0.7,
                    label=f'{strategy} (n={len(df)})', color=colors[i % len(colors)])
    
    ax1.set_xlabel('Mean Conservation Score')
    ax1.set_ylabel('Number of Genes')
    ax1.set_title('Conservation Score Distribution')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Identity distribution
    ax2 = axes[0, 1]
    if not core_df.empty:
        ax2.hist(core_df['mean_pairwise_identity'], bins=50, alpha=0.7,
                label=f'Core genes (n={len(core_df)})', color='blue')
    
    for i, (strategy, df) in enumerate(operon_dfs.items()):
        if not df.empty:
            ax2.hist(df['mean_pairwise_identity'], bins=50, alpha=0.7,
                    label=f'{strategy} (n={len(df)})', color=colors[i % len(colors)])
    
    ax2.set_xlabel('Mean Pairwise Identity')
    ax2.set_ylabel('Number of Genes')
    ax2.set_title('Sequence Identity Distribution')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Conservation vs Identity scatter
    ax3 = axes[1, 0]
    if not core_df.empty:
        ax3.scatter(core_df['mean_conservation'], core_df['mean_pairwise_identity'],
                   alpha=0.6, label='Core genes', color='blue', s=20)
    
    for i, (strategy, df) in enumerate(operon_dfs.items()):
        if not df.empty:
            ax3.scatter(df['mean_conservation'], df['mean_pairwise_identity'],
                       alpha=0.8, label=strategy, color=colors[i % len(colors)], s=50)
    
    ax3.set_xlabel('Mean Conservation Score')
    ax3.set_ylabel('Mean Pairwise Identity')
    ax3.set_title('Conservation vs Identity')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Box plot comparison
    ax4 = axes[1, 1]
    
    # Prepare data for box plot
    box_data = []
    box_labels = []
    
    if not core_df.empty:
        box_data.append(core_df['mean_conservation'])
        box_labels.append('Core genes')
    
    for strategy, df in operon_dfs.items():
        if not df.empty:
            box_data.append(df['mean_conservation'])
            box_labels.append(strategy)
    
    if box_data:
        ax4.boxplot(box_data, labels=box_labels)
        ax4.set_ylabel('Mean Conservation Score')
        ax4.set_title('Conservation Score Comparison')
        ax4.tick_params(axis='x', rotation=45)
        ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Conservation comparison plot saved to {output_file}")

def create_summary_statistics(core_df, operon_dfs, output_file):
    """Create summary statistics table."""
    stats = []
    
    # Core genes stats
    if not core_df.empty:
        stats.append({
            'dataset': 'Core genes',
            'n_genes': len(core_df),
            'mean_conservation': core_df['mean_conservation'].mean(),
            'std_conservation': core_df['mean_conservation'].std(),
            'mean_identity': core_df['mean_pairwise_identity'].mean(),
            'std_identity': core_df['mean_pairwise_identity'].std(),
            'mean_gap_percentage': core_df['gap_percentage'].mean(),
            'conservation_category': core_df['conservation_category'].mode().iloc[0] if 'conservation_category' in core_df.columns else 'Unknown'
        })
    
    # Operon strategies stats
    for strategy, df in operon_dfs.items():
        if not df.empty:
            stats.append({
                'dataset': strategy,
                'n_genes': len(df),
                'mean_conservation': df['mean_conservation'].mean(),
                'std_conservation': df['mean_conservation'].std(),
                'mean_identity': df['mean_pairwise_identity'].mean(),
                'std_identity': df['mean_pairwise_identity'].std(),
                'mean_gap_percentage': df['gap_percentage'].mean(),
                'conservation_category': df['conservation_category'].mode().iloc[0] if 'conservation_category' in df.columns else 'Unknown'
            })
    
    stats_df = pd.DataFrame(stats)
    stats_df.to_csv(output_file, index=False)
    
    print(f"✅ Summary statistics saved to {output_file}")
    
    # Print summary to console
    print("\n" + "="*60)
    print("COMPARATIVE DIVERSITY ANALYSIS SUMMARY")
    print("="*60)
    for _, row in stats_df.iterrows():
        print(f"\n{row['dataset']}:")
        print(f"  📊 Genes: {row['n_genes']}")
        print(f"  🔬 Conservation: {row['mean_conservation']:.3f} ± {row['std_conservation']:.3f}")
        print(f"  🧬 Identity: {row['mean_identity']:.3f} ± {row['std_identity']:.3f}")
        print(f"  📈 Gaps: {row['mean_gap_percentage']:.1f}%")
        print(f"  📋 Category: {row['conservation_category']}")

def main():
    parser = argparse.ArgumentParser(description="Create comparative analysis from existing CSV files")
    parser.add_argument('--core-csv', default='../04_core_gene_analysis/output/core_gene_conservation_metrics.csv',
                        help='Core gene metrics CSV')
    parser.add_argument('--strategy-a-csv', default='../05_operon_assembly_extraction/output/mappings/aa_nt_mapping/prokka/msa/operon_conservation_metrics.csv',
                        help='Strategy A operon metrics CSV')
    parser.add_argument('--strategy-d-csv', default='../05_operon_assembly_extraction/output/mappings/aa_nt_mapping/assemblies/msa/operon_conservation_metrics.csv',
                        help='Strategy D operon metrics CSV')
    parser.add_argument('--output-dir', default='output/comparative_analysis',
                        help='Output directory')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("🔍 Loading data from existing CSV files...")
    print("=" * 50)
    
    # Load data
    core_df = load_core_gene_metrics(args.core_csv)
    
    operon_dfs = {}
    if os.path.exists(args.strategy_a_csv):
        operon_dfs['Strategy A'] = load_operon_metrics(args.strategy_a_csv, 'Strategy A')
    
    if os.path.exists(args.strategy_d_csv):
        operon_dfs['Strategy D'] = load_operon_metrics(args.strategy_d_csv, 'Strategy D')
    
    if core_df.empty and not operon_dfs:
        print("❌ No data loaded. Please ensure CSV files exist.")
        return
    
    print(f"\n📊 Creating comparative analysis...")
    
    # Create comparative plots
    plot_file = os.path.join(args.output_dir, 'conservation_comparison.png')
    create_conservation_comparison_plot(core_df, operon_dfs, plot_file)
    
    # Create summary statistics
    stats_file = os.path.join(args.output_dir, 'comparative_summary.csv')
    create_summary_statistics(core_df, operon_dfs, stats_file)
    
    print(f"\n✅ Comparative analysis completed!")
    print(f"📁 Results saved to: {args.output_dir}")
    print(f"   📊 {plot_file}")
    print(f"   📋 {stats_file}")

if __name__ == "__main__":
    main()
