#!/usr/bin/env python3
"""
Plot core gene accumulation curve showing the effect of different thresholds.
Integrates with the existing core gene analysis pipeline.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
from pathlib import Path

def count_genomes_in_prokka_dir(prokka_dir):
    """Count the number of genomes in the Prokka directory."""
    genome_dirs = [d for d in Path(prokka_dir).iterdir() if d.is_dir()]
    return len(genome_dirs)

def load_gene_prevalence_data(core_genes_file, prokka_dir):
    """Load gene prevalence data from the core genes analysis."""
    
    # First, determine the total number of genomes
    n_genomes = count_genomes_in_prokka_dir(prokka_dir)
    print(f"Total genomes analyzed: {n_genomes}")
    
    # Read the core genes file to get gene counts
    if not os.path.exists(core_genes_file):
        raise FileNotFoundError(f"Core genes file not found: {core_genes_file}")
    
    core_df = pd.read_csv(core_genes_file)
    print(f"Loaded core genes data: {len(core_df)} genes")
    
    # Check if we have a 'genome_count' or 'count' column
    count_col = None
    for col in ['genome_count', 'count', 'genomes']:
        if col in core_df.columns:
            count_col = col
            break
    
    if count_col is None:
        print("Warning: No genome count column found. Will estimate from core genes file.")
        # If no count column, we'll need to create prevalence data differently
        # For now, assume all genes in the file are core (present in >= threshold genomes)
        return None, n_genomes
    
    # Create prevalence data
    prevalence_data = []
    for _, row in core_df.iterrows():
        prevalence_data.append({
            'gene': row.get('gene_name', row.get('gene', 'unknown')),
            'count': row[count_col]
        })
    
    prevalence_df = pd.DataFrame(prevalence_data)
    return prevalence_df, n_genomes

def create_threshold_analysis(prevalence_df, n_genomes, output_dir):
    """Create core gene threshold analysis plots."""
    
    if prevalence_df is None:
        print("Warning: Cannot create threshold analysis without prevalence data")
        return
    
    # Calculate core genes at different thresholds
    thresholds = np.arange(0.01, 1.01, 0.01)  # 1% to 100% in 1% steps
    core_gene_counts = []
    
    for threshold in thresholds:
        min_genomes = int(threshold * n_genomes)
        n_core = (prevalence_df['count'] >= min_genomes).sum()
        core_gene_counts.append(n_core)
        
        # Print some key thresholds
        if threshold in [0.90, 0.95, 0.99, 1.00]:
            print(f"At {threshold*100:.0f}% threshold: {n_core} core genes")
    
    # Create the main plot
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Plot the curve
    ax.plot(thresholds * 100, core_gene_counts, 'b-', linewidth=2.5, label='Core gene count')
    
    # Highlight the 95% threshold (commonly used)
    threshold_95 = 0.95
    idx_95 = int(threshold_95 * 100 - 1)
    n_core_95 = core_gene_counts[idx_95]
    
    # Add vertical line at 95%
    ax.axvline(95, color='red', linestyle='--', alpha=0.7, linewidth=2)
    ax.axhline(n_core_95, color='red', linestyle='--', alpha=0.7, linewidth=2)
    
    # Add point at 95% threshold
    ax.scatter([95], [n_core_95], color='red', s=200, zorder=5)
    ax.text(95.5, n_core_95 + max(core_gene_counts)*0.02, f'{n_core_95} genes\nat 95% threshold', 
            fontsize=12, fontweight='bold', color='red',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor='red'))
    
    # Add shaded region for "soft-core" (95-100%)
    ax.axvspan(95, 100, alpha=0.1, color='red', label='Soft-core region (≥95%)')
    
    # Add annotations for key thresholds
    for threshold, label in [(1.00, '100% (strict core)'), (0.99, '99%'), (0.90, '90%')]:
        idx = int(threshold * 100 - 1)
        n_genes = core_gene_counts[idx]
        if threshold != 0.95:  # Don't overlap with the 95% annotation
            ax.annotate(f'{n_genes} genes', 
                        xy=(threshold*100, n_genes), 
                        xytext=(threshold*100 - 5, n_genes + max(core_gene_counts)*0.05),
                        arrowprops=dict(arrowstyle='->', color='gray', alpha=0.5),
                        fontsize=10, color='gray')
    
    # Labels and title
    ax.set_xlabel('Genome prevalence threshold (%)', fontsize=14)
    ax.set_ylabel('Number of core genes', fontsize=14)
    ax.set_title(f'Core Gene Definition Across {n_genomes:,} Genomes', 
                 fontsize=16, fontweight='bold')
    
    # Grid
    ax.grid(True, alpha=0.3)
    
    # Set limits
    ax.set_xlim(0, 100)
    ax.set_ylim(0, max(core_gene_counts) * 1.1)
    
    # Add text box with key stats
    textstr = f'Total genomes: {n_genomes:,}\nTotal unique genes: {len(prevalence_df):,}'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', bbox=props)
    
    # Legend
    ax.legend(loc='upper right')
    
    # Save figure
    plt.tight_layout()
    main_plot = os.path.join(output_dir, 'core_gene_threshold_curve.png')
    plt.savefig(main_plot, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create zoomed version (85-100% range)
    if max(core_gene_counts) > 50:  # Only create zoomed view if we have enough variation
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Plot the curve (zoomed)
        ax.plot(thresholds * 100, core_gene_counts, 'b-', linewidth=2.5)
        
        # Highlight the 95% threshold
        ax.axvline(95, color='red', linestyle='--', alpha=0.7, linewidth=2)
        ax.axhline(n_core_95, color='red', linestyle='--', alpha=0.7, linewidth=2)
        ax.scatter([95], [n_core_95], color='red', s=200, zorder=5)
        
        # Add detailed annotations in zoomed view
        for pct in [85, 90, 95, 96, 97, 98, 99, 100]:
            if pct <= 100:
                idx = pct - 1
                n_genes = core_gene_counts[idx]
                ax.plot([pct, pct], [min(core_gene_counts[84:]), n_genes], 'k:', alpha=0.3)
                ax.text(pct, n_genes + (max(core_gene_counts[84:]) - min(core_gene_counts[84:]))*0.02, 
                       f'{n_genes}', ha='center', va='bottom', fontsize=9, rotation=45)
        
        # Shade the soft-core region
        ax.axvspan(95, 100, alpha=0.1, color='red', label='Soft-core region (≥95%)')
        
        # Labels
        ax.set_xlabel('Genome prevalence threshold (%)', fontsize=14)
        ax.set_ylabel('Number of core genes', fontsize=14)
        ax.set_title('Core Gene Definition (Zoomed: 85-100% Threshold)', 
                     fontsize=16, fontweight='bold')
        
        # Set limits for zoomed view
        ax.set_xlim(85, 100)
        zoom_min = min(core_gene_counts[84:])
        zoom_max = max(core_gene_counts[84:])
        margin = (zoom_max - zoom_min) * 0.1
        ax.set_ylim(zoom_min - margin, zoom_max + margin)
        
        # Grid
        ax.grid(True, alpha=0.3)
        
        # Add annotation for 95% threshold
        ax.annotate(f'{n_core_95} genes at 95%\n(commonly used threshold)', 
                    xy=(95, n_core_95), 
                    xytext=(92, n_core_95 - margin),
                    arrowprops=dict(arrowstyle='->', color='red', linewidth=2),
                    fontsize=12, fontweight='bold', color='red',
                    bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor='red'))
        
        plt.tight_layout()
        zoomed_plot = os.path.join(output_dir, 'core_gene_threshold_curve_zoomed.png')
        plt.savefig(zoomed_plot, dpi=300, bbox_inches='tight')
        plt.close()
    else:
        zoomed_plot = None
    
    # Create summary table
    summary_data = []
    for pct in [100, 99, 98, 97, 96, 95, 90, 80, 70, 50]:
        idx = pct - 1
        n_genes = core_gene_counts[idx]
        summary_data.append({
            'Threshold (%)': pct,
            'Core genes': n_genes,
            'Min genomes': int(pct/100 * n_genomes),
            'Percentage of total genes': f"{(n_genes/len(prevalence_df)*100):.1f}%"
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_file = os.path.join(output_dir, 'core_gene_threshold_summary.csv')
    summary_df.to_csv(summary_file, index=False)
    
    print(f"\nPlots and analysis saved:")
    print(f"  - {main_plot}")
    if zoomed_plot:
        print(f"  - {zoomed_plot}")
    print(f"  - {summary_file}")
    
    print(f"\nCore gene counts at key thresholds:")
    print(summary_df.to_string(index=False))
    
    return main_plot, zoomed_plot, summary_file

def main():
    parser = argparse.ArgumentParser(description="Analyze core gene thresholds")
    parser.add_argument("--core-genes-file", default="output/gene_prevalence_stats.csv",
                        help="Core genes CSV file from identify_core_genes.py")
    parser.add_argument("--prokka-dir", default="../01_prokka_annotation/output/prokka_results",
                        help="Directory containing Prokka results")
    parser.add_argument("--output-dir", default="output",
                        help="Output directory for plots and analysis")
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    try:
        # Load gene prevalence data
        prevalence_df, n_genomes = load_gene_prevalence_data(args.core_genes_file, args.prokka_dir)
        
        if prevalence_df is not None:
            # Create threshold analysis
            create_threshold_analysis(prevalence_df, n_genomes, args.output_dir)
        else:
            print("Cannot create threshold analysis without gene prevalence data.")
            print("Make sure identify_core_genes.py has been run first.")
    
    except Exception as e:
        print(f"Error: {e}")
        print("Make sure the core genes file exists and identify_core_genes.py has been run.")

if __name__ == "__main__":
    main()
