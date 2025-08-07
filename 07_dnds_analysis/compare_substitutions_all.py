#!/usr/bin/env python3
"""Compare substitution patterns between operon genes and core genes"""

import os
import glob
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import subprocess
from multiprocessing import Pool, cpu_count
import random

def get_codon_changes(codon1, codon2):
    """Classify codon change as synonymous or non-synonymous"""
    if len(codon1) != 3 or len(codon2) != 3:
        return None
    
    if '-' in codon1 or '-' in codon2 or 'N' in codon1 or 'N' in codon2:
        return None
    
    try:
        aa1 = str(Seq(codon1).translate())
        aa2 = str(Seq(codon2).translate())
        
        if aa1 == aa2:
            return 'synonymous'
        else:
            return 'non-synonymous'
    except:
        return None

def analyze_gene_substitutions(args):
    """Analyze a single gene file"""
    fasta_file, gene_type = args
    gene_name = os.path.basename(fasta_file).replace('.fasta', '').replace('_variants_aligned', '')
    
    # Files are already aligned (Steps 04 and 05), use as-is
    msa_file = fasta_file
    
    # Analyze alignment
    try:
        alignment = AlignIO.read(msa_file, 'fasta')
    except:
        return None
    
    n_seqs = len(alignment)
    if n_seqs < 2:
        return None
    
    results = {
        'gene': gene_name,
        'gene_type': gene_type,
        'n_sequences': n_seqs,
        'synonymous': 0,
        'non_synonymous': 0,
        'pairs_analyzed': 0
    }
    
    # Sample sequence pairs
    max_pairs = 20
    indices = list(range(n_seqs))
    if n_seqs > 10:
        indices = random.sample(indices, 10)
    
    # Analyze pairwise
    for i in range(len(indices)-1):
        for j in range(i+1, min(i+3, len(indices))):
            seq1 = str(alignment[indices[i]].seq).upper()
            seq2 = str(alignment[indices[j]].seq).upper()
            
            # Compare codons
            for pos in range(0, min(len(seq1), len(seq2)) - 2, 3):
                codon1 = seq1[pos:pos+3]
                codon2 = seq2[pos:pos+3]
                
                if codon1 == codon2:
                    continue
                
                change_type = get_codon_changes(codon1, codon2)
                if change_type == 'synonymous':
                    results['synonymous'] += 1
                elif change_type == 'non-synonymous':
                    results['non_synonymous'] += 1
            
            results['pairs_analyzed'] += 1
            if results['pairs_analyzed'] >= max_pairs:
                break
    
    # Calculate dN/dS
    if results['synonymous'] > 0:
        results['dn_ds_ratio'] = results['non_synonymous'] / results['synonymous']
    else:
        results['dn_ds_ratio'] = 3.0 if results['non_synonymous'] > 0 else 0
    
    # Clean up temporary alignment
    if gene_type == 'core' and os.path.exists(aligned_file):
        os.remove(aligned_file)
    
    return results

def main():
    print("=== Comparing Substitutions: Operon vs Core Genes ===\n")
    
    # Get operon files
    operon_files = glob.glob('../05_operon_assembly_extraction/output/msa/dna_alignments/*_aligned.fasta')
    print(f"Found {len(operon_files)} operon gene alignments")
    
    # Get core gene files (sample 50)
    all_core_files = glob.glob('../04_core_gene_analysis/output/core_gene_alignments/*_aligned.fasta')
    # Filter out very small files
    core_files = []
    for f in all_core_files:
        if os.path.getsize(f) > 10000:  # At least 10KB
            core_files.append(f)
    
    # Use all available core files (up to 500)
    core_files = core_files[:500]
    print(f"Using {len(core_files)} core genes for comparison")
    
    # Prepare arguments
    args_list = [(f, 'operon') for f in operon_files] + [(f, 'core') for f in core_files]
    
    # Process in parallel
    print("\nAnalyzing genes...")
    n_cores = min(cpu_count() - 1, 8)
    with Pool(n_cores) as pool:
        results = pool.map(analyze_gene_substitutions, args_list)
    
    # Filter out None results
    results = [r for r in results if r is not None]
    df = pd.DataFrame(results)
    
    print(f"\nSuccessfully analyzed {len(df)} genes")
    print(f"  Operon genes: {len(df[df['gene_type'] == 'operon'])}")
    print(f"  Core genes: {len(df[df['gene_type'] == 'core'])}")
    
    # Create comparison plots
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Plot 1: dN/dS distribution comparison
    ax1 = axes[0, 0]
    operon_dnds = df[df['gene_type'] == 'operon']['dn_ds_ratio'].values
    core_dnds = df[df['gene_type'] == 'core']['dn_ds_ratio'].values
    
    # Cap extreme values for visualization
    operon_dnds = np.clip(operon_dnds, 0, 3)
    core_dnds = np.clip(core_dnds, 0, 3)
    
    bp = ax1.boxplot([operon_dnds, core_dnds], positions=[1, 2], widths=0.6, 
                     patch_artist=True, showmeans=True)
    bp['boxes'][0].set_facecolor('red')
    bp['boxes'][1].set_facecolor('blue')
    
    ax1.set_xticklabels(['Operon genes', 'Core genes'])
    ax1.set_ylabel('dN/dS ratio', fontsize=12)
    ax1.set_title('Selection Pressure Comparison', fontsize=14, fontweight='bold')
    ax1.axhline(y=1, color='gray', linestyle='--', alpha=0.5)
    ax1.text(1.5, 1.05, 'Neutral evolution', ha='center', fontsize=10, color='gray')
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Add statistical test
    from scipy import stats
    if len(operon_dnds) > 0 and len(core_dnds) > 0:
        statistic, pvalue = stats.mannwhitneyu(operon_dnds, core_dnds, alternative='two-sided')
        ax1.text(0.02, 0.95, f'Mann-Whitney U test\np-value: {pvalue:.3e}', 
                transform=ax1.transAxes, fontsize=10,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Plot 2: Histogram of dN/dS values
    ax2 = axes[0, 1]
    bins = np.linspace(0, 2, 21)
    ax2.hist(operon_dnds[operon_dnds < 2], bins=bins, alpha=0.5, label='Operon', color='red', density=True)
    ax2.hist(core_dnds[core_dnds < 2], bins=bins, alpha=0.5, label='Core', color='blue', density=True)
    ax2.axvline(x=1, color='gray', linestyle='--', alpha=0.5)
    ax2.set_xlabel('dN/dS ratio', fontsize=12)
    ax2.set_ylabel('Density', fontsize=12)
    ax2.set_title('Distribution of dN/dS Values', fontsize=14, fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Scatter plot
    ax3 = axes[1, 0]
    operon_df = df[df['gene_type'] == 'operon']
    core_df = df[df['gene_type'] == 'core']
    
    # Plot core genes first (blue, underneath)
    if len(core_df) > 0:
        ax3.scatter(core_df['synonymous'], core_df['non_synonymous'], 
                   color='blue', alpha=0.6, s=50, label='Core genes', zorder=1)
    # Then plot operon genes on top (red, above)
    if len(operon_df) > 0:
        ax3.scatter(operon_df['synonymous'], operon_df['non_synonymous'], 
                   color='red', alpha=0.8, s=50, label='Operon genes', zorder=2, edgecolors='darkred', linewidth=1)
    
    # Add diagonal lines for different dN/dS ratios
    max_val = max(df['synonymous'].max(), df['non_synonymous'].max())
    x_line = np.linspace(0, max_val, 100)
    ax3.plot(x_line, x_line, 'k--', alpha=0.3, label='dN/dS = 1')
    ax3.plot(x_line, 0.5*x_line, 'g--', alpha=0.3, label='dN/dS = 0.5')
    ax3.plot(x_line, 0.1*x_line, 'b--', alpha=0.3, label='dN/dS = 0.1')
    
    ax3.set_xlabel('Synonymous substitutions', fontsize=12)
    ax3.set_ylabel('Non-synonymous substitutions', fontsize=12)
    ax3.set_title('Substitution Counts', fontsize=14, fontweight='bold')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Summary statistics
    ax4 = axes[1, 1]
    ax4.axis('off')
    
    # Calculate statistics
    operon_stats = df[df['gene_type'] == 'operon']['dn_ds_ratio'].describe()
    core_stats = df[df['gene_type'] == 'core']['dn_ds_ratio'].describe()
    
    summary_text = "Summary Statistics\n\n"
    summary_text += "Operon genes:\n"
    summary_text += f"  • n = {len(operon_df)}\n"
    summary_text += f"  • Mean dN/dS: {operon_stats['mean']:.3f}\n"
    summary_text += f"  • Median dN/dS: {operon_stats['50%']:.3f}\n"
    summary_text += f"  • Range: {operon_stats['min']:.3f} - {operon_stats['max']:.3f}\n"
    
    summary_text += "\nCore genes:\n"
    summary_text += f"  • n = {len(core_df)}\n"
    summary_text += f"  • Mean dN/dS: {core_stats['mean']:.3f}\n"
    summary_text += f"  • Median dN/dS: {core_stats['50%']:.3f}\n"
    summary_text += f"  • Range: {core_stats['min']:.3f} - {core_stats['max']:.3f}\n"
    
    # Count genes under strong selection
    strong_selection_operon = (operon_df['dn_ds_ratio'] < 0.5).sum()
    strong_selection_core = (core_df['dn_ds_ratio'] < 0.5).sum()
    
    summary_text += f"\nStrong purifying selection (dN/dS < 0.5):\n"
    summary_text += f"  • Operon: {strong_selection_operon}/{len(operon_df)} ({strong_selection_operon/len(operon_df)*100:.0f}%)\n"
    summary_text += f"  • Core: {strong_selection_core}/{len(core_df)} ({strong_selection_core/len(core_df)*100:.0f}%)\n"
    
    ax4.text(0.1, 0.9, summary_text, transform=ax4.transAxes, fontsize=12,
             verticalalignment='top', 
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.suptitle('Operon vs Core Genes: Selection Pressure Analysis', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig('operon_vs_core_substitutions.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save detailed results
    df.to_csv('substitution_comparison_results.csv', index=False)
    
    print("\n=== Results Summary ===")
    print(f"Operon genes mean dN/dS: {operon_stats['mean']:.3f}")
    print(f"Core genes mean dN/dS: {core_stats['mean']:.3f}")
    print(f"\nFiles saved:")
    print("  - operon_vs_core_substitutions.png")
    print("  - substitution_comparison_results.csv")

if __name__ == '__main__':
    main()