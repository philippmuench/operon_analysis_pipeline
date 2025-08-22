#!/usr/bin/env python3
"""
Comparative analysis of core genes vs operon genes conservation.
Creates bar plot visualizations showing conservation rankings.
Also creates conservation plots from MSA alignments.
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from typing import Dict, List, Tuple, Optional
import argparse

# Set style for publication-quality figures
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

def load_core_gene_metrics():
    """Load pre-computed core gene conservation metrics from step 04."""
    core_csv = "../04_core_gene_analysis/output/core_gene_conservation_metrics.csv"
    
    if not os.path.exists(core_csv):
        print(f"ERROR: Core gene metrics not found at {core_csv}")
        print("   Please run step 04 first: cd ../04_core_gene_analysis && sbatch run_core_analysis.sh")
        return None
    
    df = pd.read_csv(core_csv)
    print(f"Loaded {len(df)} core genes from step 04")
    
    # Fix corrupted column name in core gene CSV (appears as 'in' instead of 'mean_conservation')
    if 'in' in df.columns and 'mean_conservation' not in df.columns:
        print("   NOTE: Fixing corrupted column name 'in' -> 'mean_conservation'")
        df.rename(columns={'in': 'mean_conservation'}, inplace=True)
    
    # Load gene descriptions from prevalence stats if available
    prevalence_csv = "../04_core_gene_analysis/output/gene_prevalence_stats.csv"
    if os.path.exists(prevalence_csv):
        prevalence_df = pd.read_csv(prevalence_csv)
        # Create a mapping of gene to product description
        gene_to_product = dict(zip(prevalence_df['gene'], prevalence_df['product']))
        
        # Add gene descriptions to the dataframe
        df['gene_full_name'] = df['gene'].apply(
            lambda g: f"{g} ({gene_to_product[g]})" if g in gene_to_product and pd.notna(gene_to_product[g])
            else g
        )
        print("   Loaded gene descriptions from prevalence stats")
    else:
        # Use gene_full_name if available, otherwise construct it
        if 'gene_full_name' not in df.columns:
            if 'gene_description' in df.columns:
                df['gene_full_name'] = df.apply(
                    lambda row: f"{row['gene']} ({row['gene_description']})" 
                    if pd.notna(row.get('gene_description')) and row['gene'] != row['gene_description']
                    else row['gene'],
                    axis=1
                )
            else:
                df['gene_full_name'] = df['gene']
    
    # Ensure required columns exist
    required_cols = ['gene', 'mean_conservation']
    if not all(col in df.columns for col in required_cols):
        print(f"   Available columns: {list(df.columns)}")
        print(f"   WARNING: Missing required columns: {[col for col in required_cols if col not in df.columns]}")
    
    return df

def load_operon_metrics():
    """Load pre-computed operon conservation metrics from step 05."""
    operon_metrics = {}
    
    # Define paths to different extraction strategies - use streamlined output
    # Try detailed metrics first (matches core gene format), fallback to basic if not available
    detailed_path = '../05_operon_assembly_extraction/output/operon_conservation_metrics_detailed.csv'
    basic_path = '../05_operon_assembly_extraction/output/operon_conservation_metrics.csv'
    
    strategies = {
        'Operon Genes': detailed_path if os.path.exists(detailed_path) else basic_path,
    }
    
    # Also check for alternative strategies if needed (legacy paths)
    alternative_strategies = {
        'Primary (Assembly-tblastn)': '../05_operon_assembly_extraction/output/mappings/aa_nt_mapping/assemblies/msa/operon_conservation_metrics.csv',
        'Prokka-tblastn': '../05_operon_assembly_extraction/output/mappings/aa_nt_mapping/prokka/msa/operon_conservation_metrics.csv',
        'Prokka-blastn': '../05_operon_assembly_extraction/output/mappings/nt_nt_mapping/prokka_genome/msa/operon_conservation_metrics.csv'
    }
    
    # Load primary strategy
    for strategy_name, csv_path in strategies.items():
        if os.path.exists(csv_path):
            df = pd.read_csv(csv_path)
            operon_metrics[strategy_name] = df
            print(f"Loaded {len(df)} operon genes from {strategy_name}")
        else:
            print(f"ERROR: Primary operon metrics not found at {csv_path}")
            print("   Checking alternative strategies...")
            
            # Try alternative strategies
            for alt_name, alt_path in alternative_strategies.items():
                if os.path.exists(alt_path):
                    df = pd.read_csv(alt_path)
                    operon_metrics[alt_name] = df
                    print(f"Loaded {len(df)} operon genes from {alt_name}")
                    break
    
    if not operon_metrics:
        print("ERROR: No operon metrics found. Please run step 05 first.")
        return None
    
    # Standardize column names for compatibility
    for strategy_name, df in operon_metrics.items():
        # Rename columns to match expected names
        column_mapping = {
            'pairwise_identity': 'mean_pairwise_identity',
            'conservation_score': 'mean_conservation',
            'num_sequences': 'n_sequences',
            'gap_percentage': 'mean_gaps_per_position'
        }
        df.rename(columns=column_mapping, inplace=True)
        
        # Add gene_full_name if not present
        if 'gene_full_name' not in df.columns:
            df['gene_full_name'] = df['gene']
    
    return operon_metrics

def create_conservation_bar_plot(core_df, operon_df, output_dir):
    """Create bar plot showing top, operon, and bottom conserved genes."""
    
    # Use mean_pairwise_identity if available, otherwise use mean_conservation
    if 'mean_pairwise_identity' in core_df.columns:
        metric_col = 'mean_pairwise_identity'
        metric_name = 'Average Pairwise Identity (%)'
        # Convert to percentage if needed
        if core_df[metric_col].max() <= 1:
            core_df[metric_col] = core_df[metric_col] * 100
            if operon_df is not None and metric_col in operon_df.columns:
                operon_df[metric_col] = operon_df[metric_col] * 100
    else:
        metric_col = 'mean_conservation'
        metric_name = 'Conservation Score (%)'
        # Convert to percentage for mean_conservation too (0-1 to 0-100)
        if core_df[metric_col].max() <= 1:
            core_df[metric_col] = core_df[metric_col] * 100
            if operon_df is not None and metric_col in operon_df.columns:
                operon_df[metric_col] = operon_df[metric_col] * 100
    
    # Sort core genes by conservation
    core_df_sorted = core_df.sort_values(metric_col, ascending=False)
    
    # Get operon gene data and calculate ranks
    operon_data = []
    if operon_df is not None and metric_col in operon_df.columns:
        for idx, row in operon_df.iterrows():
            gene_name = row['gene']
            conservation_value = row[metric_col]
            
            # Calculate rank among all core genes
            rank = (core_df[metric_col] > conservation_value).sum() + 1
            percentile = (len(core_df) - rank + 1) / len(core_df) * 100
            
            operon_data.append({
                'gene': gene_name,
                'conservation': conservation_value,
                'rank': rank,
                'percentile': percentile
            })
    
    # Sort operon genes by conservation
    operon_sorted = sorted(operon_data, key=lambda x: x['conservation'], reverse=True)
    
    # Create figure with three panels - wider to accommodate longer gene names
    fig = plt.figure(figsize=(16, 16), facecolor='white')
    gs = fig.add_gridspec(3, 1, height_ratios=[20, 7, 20], hspace=0.3)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    ax3 = fig.add_subplot(gs[2])
    
    # Set background color
    for ax in [ax1, ax2, ax3]:
        ax.set_facecolor('#f8f8f8')
    
    # Determine x-axis limits based on data range
    all_values = core_df[metric_col].values
    if operon_data:
        all_values = np.concatenate([all_values, [d['conservation'] for d in operon_data]])
    x_min = max(0, np.min(all_values) - 5)
    x_max = min(100, np.max(all_values) + 2)
    
    # Top 20 most conserved genes
    top20 = core_df_sorted.head(20).iloc[::-1]  # Reverse for display
    bars1 = ax1.barh(range(20), top20[metric_col], color='darkgreen', alpha=0.7, edgecolor='black', linewidth=0.5)
    
    # Add labels
    for i, (bar, (idx, row)) in enumerate(zip(bars1, top20.iterrows())):
        rank = (core_df[metric_col] > row[metric_col]).sum() + 1
        percentile = (len(core_df) - rank + 1) / len(core_df) * 100
        # Add text inside bar
        ax1.text(bar.get_width() - 1, bar.get_y() + bar.get_height()/2, 
                 f'{row[metric_col]:.1f}% | #{rank}', 
                 va='center', ha='right', fontsize=8, color='white', fontweight='bold')
    
    ax1.set_yticks(range(20))
    # Use full gene names if available, otherwise use gene names
    top20_labels = []
    for idx, row in top20.iterrows():
        if 'gene_full_name' in row and pd.notna(row['gene_full_name']):
            gene_display = str(row['gene_full_name'])
        else:
            gene_display = str(row['gene']).replace('_', ' ')
        # Allow up to 80 characters for better readability
        top20_labels.append(gene_display[:80] + '...' if len(gene_display) > 80 else gene_display)
    ax1.set_yticklabels(top20_labels, fontsize=8)
    ax1.set_xlabel(metric_name, fontsize=11)
    ax1.set_title('Top 20 Most Conserved Core Genes', fontsize=13, fontweight='bold')
    ax1.set_xlim(x_min, x_max)
    ax1.grid(axis='x', alpha=0.3)
    
    # Operon genes panel
    if operon_sorted:
        bars2 = ax2.barh(range(len(operon_sorted)), 
                         [d['conservation'] for d in operon_sorted], 
                         color='red', alpha=0.8, edgecolor='black', linewidth=0.5)
        
        # Add labels
        for i, (bar, d) in enumerate(zip(bars2, operon_sorted)):
            ax2.text(bar.get_width() - 1, bar.get_y() + bar.get_height()/2, 
                     f'{d["conservation"]:.1f}% | #{d["rank"]} ({d["percentile"]:.0f}%ile)', 
                     va='center', ha='right', fontsize=8, color='white', fontweight='bold')
        
        ax2.set_yticks(range(len(operon_sorted)))
        # Use full descriptive gene names for operon genes
        operon_labels = []
        gene_name_mapping = {
            'frpC': 'Fructoselysine-6-phosphate deglycase (frpC)',
            'glpC': 'Glucoselysine-6-phosphate deglycase (glpC)',
            'ptsD': 'PTS system fructoselysine/glucoselysine IID component (ptsD)',
            'ptsC': 'PTS system fructoselysine/glucoselysine IIC component (ptsC)',
            'ptsB': 'PTS system fructoselysine/glucoselysine IIB component (ptsB)',
            'ptsA': 'PTS system fructoselysine/glucoselysine IIA component (ptsA)',
            'fruR': 'Sigma-54 dependent transcriptional regulator (fruR)'
        }
        
        for d in operon_sorted:
            gene_name = d['gene']
            # Check if this is a known gene name and map it
            if gene_name in gene_name_mapping:
                display_name = gene_name_mapping[gene_name]
            else:
                # Otherwise, clean up the gene name
                display_name = gene_name.replace('_', ' ')
                # Handle common patterns
                if 'fructoselysine' in display_name.lower() or 'glucoselysine' in display_name.lower():
                    display_name = display_name.replace('fructoselysine glucoselysine PTS system', 'PTS system fructoselysine/glucoselysine')
                if 'sigma-54' in display_name.lower():
                    display_name = display_name.replace('sigma-54 dependent transcriptional regulator', 'Sigma-54 dependent transcriptional regulator')
            
            operon_labels.append(display_name)
        
        ax2.set_yticklabels(operon_labels, fontsize=7.5)
        ax2.set_xlabel(metric_name, fontsize=11)
        ax2.set_title('Fructoselysine/Glucoselysine Operon Genes', fontsize=13, fontweight='bold', color='darkred')
        ax2.set_xlim(x_min, x_max)
        ax2.grid(axis='x', alpha=0.3)
        ax2.invert_yaxis()
    
    # Bottom 20 least conserved genes
    bottom20 = core_df_sorted.tail(20).iloc[::-1]  # Reverse for display
    bars3 = ax3.barh(range(20), bottom20[metric_col], color='darkred', alpha=0.7, edgecolor='black', linewidth=0.5)
    
    # Add labels
    for i, (bar, (idx, row)) in enumerate(zip(bars3, bottom20.iterrows())):
        rank = (core_df[metric_col] > row[metric_col]).sum() + 1
        percentile = (len(core_df) - rank + 1) / len(core_df) * 100
        # Add text inside bar
        if bar.get_width() > 10:  # Only add text if bar is wide enough
            ax3.text(bar.get_width() - 1, bar.get_y() + bar.get_height()/2, 
                     f'{row[metric_col]:.1f}% | #{rank}', 
                     va='center', ha='right', fontsize=8, color='white', fontweight='bold')
        else:
            ax3.text(bar.get_width() + 1, bar.get_y() + bar.get_height()/2, 
                     f'{row[metric_col]:.1f}% | #{rank}', 
                     va='center', ha='left', fontsize=8, color='black')
    
    ax3.set_yticks(range(20))
    # Use full gene names if available, otherwise use gene names
    bottom20_labels = []
    for idx, row in bottom20.iterrows():
        if 'gene_full_name' in row and pd.notna(row['gene_full_name']):
            gene_display = str(row['gene_full_name'])
        else:
            gene_display = str(row['gene']).replace('_', ' ')
        # Allow up to 80 characters for better readability
        bottom20_labels.append(gene_display[:80] + '...' if len(gene_display) > 80 else gene_display)
    ax3.set_yticklabels(bottom20_labels, fontsize=8)
    ax3.set_xlabel(metric_name, fontsize=11)
    ax3.set_title('Bottom 20 Least Conserved Core Genes', fontsize=13, fontweight='bold')
    ax3.set_xlim(x_min, x_max)
    ax3.grid(axis='x', alpha=0.3)
    
    # Add reference lines
    for ax in [ax1, ax2, ax3]:
        for val in [50, 75, 90, 95, 99]:
            if x_min <= val <= x_max:
                ax.axvline(val, color='gray', linestyle=':', alpha=0.3, linewidth=0.5)
    
    plt.suptitle('Core Gene Conservation Analysis', fontsize=15, fontweight='bold', y=0.995)
    plt.tight_layout()
    
    # Save figure
    output_file = os.path.join(output_dir, 'conservation_ranking_bars.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(os.path.join(output_dir, 'conservation_ranking_bars.pdf'), bbox_inches='tight', facecolor='white')
    print(f"Saved conservation ranking plot to {output_file}")
    
    plt.close()
    
    return operon_data

def create_distribution_plot(core_df, operon_data, output_dir):
    """Create histogram showing distribution of conservation with operon genes highlighted."""
    
    fig, ax = plt.subplots(figsize=(12, 7), facecolor='white')
    
    # Use appropriate metric column
    if 'mean_pairwise_identity' in core_df.columns:
        metric_col = 'mean_pairwise_identity'
        metric_name = 'Average Pairwise Identity (%)'
    else:
        metric_col = 'mean_conservation'
        metric_name = 'Conservation Score'
    
    # Create histogram
    n, bins, patches = ax.hist(core_df[metric_col], bins=50, alpha=0.6, 
                                color='gray', edgecolor='black', linewidth=0.5)
    
    # Highlight operon genes if available
    if operon_data:
        for i, gene_data in enumerate(operon_data):
            ax.axvline(gene_data['conservation'], color='red', linestyle='--', 
                      alpha=0.8, linewidth=2)
            
            # Add label
            y_pos = ax.get_ylim()[1] * (0.9 - i * 0.08)
            gene_short = gene_data['gene'].replace('_', ' ')
            if len(gene_short) > 30:
                gene_short = gene_short[:27] + '...'
            
            ax.text(gene_data['conservation'] + 0.5, y_pos, 
                   f"{gene_short}\nRank: {gene_data['rank']}/{len(core_df)}",
                   fontsize=8, ha='left', 
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7))
    
    ax.set_xlabel(metric_name, fontsize=12)
    ax.set_ylabel('Number of Core Genes', fontsize=12)
    ax.set_title('Distribution of Core Gene Conservation with Operon Genes Highlighted', 
                fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'conservation_distribution.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(os.path.join(output_dir, 'conservation_distribution.pdf'), bbox_inches='tight', facecolor='white')
    print(f"Saved distribution plot to {output_file}")
    
    plt.close()

def save_summary_table(core_df, operon_data, output_dir):
    """Save summary statistics to CSV file."""
    
    if not operon_data:
        print("No operon data available for summary table")
        return
    
    # Create summary dataframe
    summary_data = []
    for gene_data in operon_data:
        summary_data.append({
            'Gene': gene_data['gene'],
            'Conservation (%)': f"{gene_data['conservation']:.2f}",
            'Rank': gene_data['rank'],
            'Total Core Genes': len(core_df),
            'Percentile': f"{gene_data['percentile']:.1f}"
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df = summary_df.sort_values('Rank')
    
    # Save to CSV
    output_file = os.path.join(output_dir, 'operon_conservation_summary.csv')
    summary_df.to_csv(output_file, index=False)
    print(f"Saved summary table to {output_file}")
    
    # Print summary statistics
    print("\n" + "="*60)
    print("CONSERVATION ANALYSIS SUMMARY")
    print("="*60)
    
    # Use appropriate metric
    if 'mean_pairwise_identity' in core_df.columns:
        metric_col = 'mean_pairwise_identity'
        metric_name = 'pairwise identity'
    else:
        metric_col = 'mean_conservation'
        metric_name = 'conservation'
    
    print(f"\nCore Genes (n={len(core_df)}):")
    print(f"  Mean {metric_name}: {core_df[metric_col].mean():.2f}%")
    print(f"  Median {metric_name}: {core_df[metric_col].median():.2f}%")
    print(f"  Range: {core_df[metric_col].min():.2f}% - {core_df[metric_col].max():.2f}%")
    
    if operon_data:
        operon_conservations = [d['conservation'] for d in operon_data]
        print(f"\nOperon Genes (n={len(operon_data)}):")
        print(f"  Mean {metric_name}: {np.mean(operon_conservations):.2f}%")
        print(f"  Median {metric_name}: {np.median(operon_conservations):.2f}%")
        print(f"  Range: {np.min(operon_conservations):.2f}% - {np.max(operon_conservations):.2f}%")
        
        # Compare to core genes
        difference = np.mean(operon_conservations) - core_df[metric_col].mean()
        if difference > 0:
            print(f"\nOperon genes are {difference:.2f}% MORE conserved than average core genes")
        else:
            print(f"\nOperon genes are {abs(difference):.2f}% LESS conserved than average core genes")
        
        # Show individual operon gene rankings
        print(f"\nOperon Gene Rankings:")
        for gene_data in sorted(operon_data, key=lambda x: x['rank']):
            print(f"  {gene_data['gene'][:40]:40s} : Rank {gene_data['rank']:4d}/{len(core_df)} ({gene_data['percentile']:.0f} percentile)")
    
    print("="*60)

def calculate_shannon_entropy(alignment_file: str) -> Tuple[List[float], Dict]:
    """
    Calculate Shannon entropy (conservation score) for each position in an alignment.
    
    Args:
        alignment_file: Path to alignment file
    
    Returns:
        Tuple of conservation scores and metadata
    """
    try:
        alignment = AlignIO.read(alignment_file, "fasta")
    except Exception as e:
        print(f"    Error reading alignment: {e}")
        return [], {}
    
    if len(alignment) == 0:
        return [], {}
    
    alignment_length = alignment.get_alignment_length()
    conservation_scores = []
    
    for pos in range(alignment_length):
        # Get column
        column = alignment[:, pos]
        
        # Count frequencies (excluding gaps)
        freq_dict = {}
        total = 0
        for base in column:
            if base not in ['-', 'N', 'n']:
                freq_dict[base] = freq_dict.get(base, 0) + 1
                total += 1
        
        # Calculate Shannon entropy
        if total > 0:
            entropy = 0
            for count in freq_dict.values():
                freq = count / total
                if freq > 0:
                    entropy -= freq * np.log2(freq)
            
            # Normalize by max possible entropy (log2(4) for DNA)
            max_entropy = np.log2(min(4, len(freq_dict))) if len(freq_dict) > 1 else 0
            conservation = 1 - (entropy / max_entropy) if max_entropy > 0 else 1
        else:
            conservation = 0  # All gaps
        
        conservation_scores.append(conservation)
    
    metadata = {
        'num_sequences': len(alignment),
        'alignment_length': alignment_length,
        'mean_conservation': np.mean(conservation_scores) if conservation_scores else 0
    }
    
    return conservation_scores, metadata

def create_conservation_plots(msa_dir: str, output_dir: str, 
                            title_suffix: str = "") -> None:
    """
    Create enhanced conservation plots with Shannon entropy.
    
    Args:
        msa_dir: Directory containing alignment files
        output_dir: Directory for output plots
        title_suffix: Suffix for plot titles
    """
    print(f"\n{'='*60}")
    print("Creating conservation plots")
    print(f"{'='*60}")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Find alignment files
    align_files = [f for f in os.listdir(msa_dir) 
                   if f.endswith('_aligned.fasta') or f.endswith('.fasta')]
    
    if not align_files:
        print(f"Warning: No alignment files found in {msa_dir}")
        return
    
    for align_file in align_files:
        align_path = os.path.join(msa_dir, align_file)
        gene_name = align_file.replace('_aligned.fasta', '').replace('.fasta', '')
        
        print(f"  Processing {gene_name}...")
        
        # Calculate conservation
        conservation_scores, metadata = calculate_shannon_entropy(align_path)
        
        if not conservation_scores:
            print(f"    Warning: No valid alignment for {gene_name}")
            continue
        
        # Create Shannon entropy plot
        fig, ax = plt.subplots(figsize=(15, 6))
        positions = np.arange(len(conservation_scores))
        
        ax.plot(positions, conservation_scores, color='blue', linewidth=1.5, alpha=0.8)
        ax.fill_between(positions, conservation_scores, alpha=0.3, color='blue')
        
        ax.set_xlabel('Position', fontsize=12)
        ax.set_ylabel('Conservation Score (Shannon Entropy)', fontsize=12)
        ax.set_title(f'{gene_name} Conservation{" - " + title_suffix if title_suffix else ""}',
                    fontsize=14, fontweight='bold')
        ax.set_ylim(0, 1.05)
        ax.grid(True, alpha=0.3)
        
        # Add metadata to plot
        info_text = f"Sequences: {metadata['num_sequences']}\nLength: {metadata['alignment_length']} bp\nMean Conservation: {metadata['mean_conservation']:.3f}"
        ax.text(0.02, 0.98, info_text, transform=ax.transAxes, 
               fontsize=10, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        # Add conservation threshold lines
        ax.axhline(y=0.9, color='green', linestyle='--', alpha=0.5, label='High (>0.9)')
        ax.axhline(y=0.7, color='orange', linestyle='--', alpha=0.5, label='Medium (>0.7)')
        ax.axhline(y=0.5, color='red', linestyle='--', alpha=0.5, label='Low (>0.5)')
        ax.legend(loc='upper right', fontsize=10)
        
        plt.tight_layout()
        
        # Save plot
        output_file = os.path.join(output_dir, f"{gene_name}_conservation.png")
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.savefig(output_file.replace('.png', '.pdf'), bbox_inches='tight')
        plt.close()
        
        print(f"    Saved: {output_file}")

def main():
    """Main analysis pipeline."""
    
    parser = argparse.ArgumentParser(description='Comparative conservation analysis and plotting')
    parser.add_argument('--mode', choices=['comparative', 'plots', 'both'], default='comparative',
                       help='Analysis mode: comparative (bar plots), plots (conservation plots), or both')
    parser.add_argument('--msa-dir', help='Directory with MSA files (for plots mode)')
    parser.add_argument('--plot-dir', help='Output directory for conservation plots')
    parser.add_argument('--title-suffix', default='', help='Suffix for plot titles')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)
    
    if args.mode in ['comparative', 'both']:
        print("="*60)
        print("COMPARATIVE CONSERVATION ANALYSIS")
        print("Core Genes vs Operon Genes")
        print("="*60 + "\n")
        
        # Load pre-computed metrics
        core_df = load_core_gene_metrics()
        operon_dfs = load_operon_metrics()
        
        if not operon_dfs:
            print("\nERROR: Cannot proceed without operon metrics. Please run step 05 first.")
            if args.mode == 'comparative':
                sys.exit(1)
        else:
            # Use the first (primary) strategy for analysis
            primary_strategy = list(operon_dfs.keys())[0]
            operon_df = operon_dfs[primary_strategy]
            print(f"\nUsing {primary_strategy} for analysis")
            
            print("\n" + "="*60)
            print("GENERATING COMPARATIVE VISUALIZATIONS...")
            print("="*60 + "\n")
            
            # Create conservation ranking bar plot
            operon_data = create_conservation_bar_plot(core_df, operon_df, output_dir)
            
            # Create distribution plot
            create_distribution_plot(core_df, operon_data, output_dir)
            
            # Save summary table and statistics
            save_summary_table(core_df, operon_data, output_dir)
    
    if args.mode in ['plots', 'both']:
        # Generate conservation plots from MSA
        if args.msa_dir:
            plot_dir = args.plot_dir or os.path.join(output_dir, 'conservation_plots')
            create_conservation_plots(args.msa_dir, plot_dir, args.title_suffix)
        else:
            # Default to looking for MSA files from step 05
            default_msa = "../05_operon_assembly_extraction/output/msa/dna_alignments"
            if os.path.exists(default_msa):
                plot_dir = args.plot_dir or os.path.join(output_dir, 'conservation_plots')
                print(f"\nUsing default MSA directory: {default_msa}")
                create_conservation_plots(default_msa, plot_dir, args.title_suffix)
            else:
                print("\nERROR: No MSA directory specified or found. Use --msa-dir option.")
    
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE!")
    print(f"All results saved to: {output_dir}/")
    if args.mode in ['comparative', 'both']:
        print("Generated comparative files:")
        print("  - conservation_ranking_bars.png/pdf")
        print("  - conservation_distribution.png/pdf")
        print("  - operon_conservation_summary.csv")
    if args.mode in ['plots', 'both']:
        print("Generated conservation plots:")
        print("  - conservation_plots/*.png")
        print("  - conservation_plots/*.pdf")
    print("="*60)

if __name__ == "__main__":
    main()