#!/usr/bin/env python3
"""
Compare dN/dS results between operon and core genes.
"""

import pandas as pd
import numpy as np
from scipy import stats
import argparse
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--operon_file', required=True)
    parser.add_argument('--core_file', required=True)
    parser.add_argument('--output', required=True)
    
    args = parser.parse_args()
    
    # Load data
    operon_df = pd.read_csv(args.operon_file)
    core_df = pd.read_csv(args.core_file)
    
    # Filter out infinite values
    operon_dnds = operon_df['dn_ds_ratio'][operon_df['dn_ds_ratio'] != float('inf')]
    core_dnds = core_df['dn_ds_ratio'][core_df['dn_ds_ratio'] != float('inf')]
    
    # Calculate statistics
    operon_stats = {
        'n': len(operon_dnds),
        'mean': operon_dnds.mean(),
        'median': operon_dnds.median(),
        'std': operon_dnds.std(),
        'min': operon_dnds.min(),
        'max': operon_dnds.max()
    }
    
    core_stats = {
        'n': len(core_dnds),
        'mean': core_dnds.mean(),
        'median': core_dnds.median(),
        'std': core_dnds.std(),
        'min': core_dnds.min(),
        'max': core_dnds.max()
    }
    
    # Statistical test
    statistic, pvalue = stats.mannwhitneyu(operon_dnds, core_dnds, alternative='less')
    
    # Selection categories
    operon_selection = {
        'strong_purifying': (operon_dnds < 0.5).sum(),
        'weak_purifying': ((operon_dnds >= 0.5) & (operon_dnds < 1)).sum(),
        'neutral': ((operon_dnds >= 0.9) & (operon_dnds <= 1.1)).sum(),
        'positive': (operon_dnds > 1.5).sum()
    }
    
    core_selection = {
        'strong_purifying': (core_dnds < 0.5).sum(),
        'weak_purifying': ((core_dnds >= 0.5) & (core_dnds < 1)).sum(),
        'neutral': ((core_dnds >= 0.9) & (core_dnds <= 1.1)).sum(),
        'positive': (core_dnds > 1.5).sum()
    }
    
    # Write summary
    with open(args.output, 'w') as f:
        f.write("=== dN/dS Comparison: Operon vs Core Genes ===\n\n")
        
        f.write("Operon Genes:\n")
        f.write(f"  N genes: {operon_stats['n']}\n")
        f.write(f"  Mean dN/dS: {operon_stats['mean']:.3f}\n")
        f.write(f"  Median dN/dS: {operon_stats['median']:.3f}\n")
        f.write(f"  Range: {operon_stats['min']:.3f} - {operon_stats['max']:.3f}\n")
        f.write("\n  Selection categories:\n")
        for cat, count in operon_selection.items():
            f.write(f"    {cat}: {count} ({count/operon_stats['n']*100:.1f}%)\n")
        
        f.write("\n\nCore Genes:\n")
        f.write(f"  N genes: {core_stats['n']}\n")
        f.write(f"  Mean dN/dS: {core_stats['mean']:.3f}\n")
        f.write(f"  Median dN/dS: {core_stats['median']:.3f}\n")
        f.write(f"  Range: {core_stats['min']:.3f} - {core_stats['max']:.3f}\n")
        f.write("\n  Selection categories:\n")
        for cat, count in core_selection.items():
            f.write(f"    {cat}: {count} ({count/core_stats['n']*100:.1f}%)\n")
        
        f.write(f"\n\nStatistical Comparison:\n")
        f.write(f"  Mann-Whitney U test (operon < core): p = {pvalue:.3e}\n")
        
        if pvalue < 0.05:
            f.write("  Result: Operon genes have significantly lower dN/dS than core genes\n")
        else:
            f.write("  Result: No significant difference in dN/dS\n")
        
        f.write("\n\nInterpretation:\n")
        if operon_stats['mean'] < 0.5:
            f.write("  - Operon genes are under STRONG purifying selection\n")
        elif operon_stats['mean'] < 1:
            f.write("  - Operon genes are under purifying selection\n")
        else:
            f.write("  - Operon genes show neutral or positive selection\n")
        
        f.write(f"  - Operon genes are {core_stats['mean']/operon_stats['mean']:.1f}x more conserved than average core genes\n")
    
    # Create visualization
    plt.figure(figsize=(10, 6))
    
    # Box plot
    plt.subplot(1, 2, 1)
    data = [operon_dnds, core_dnds]
    labels = ['Operon\nGenes', 'Core\nGenes']
    bp = plt.boxplot(data, labels=labels, patch_artist=True)
    bp['boxes'][0].set_facecolor('red')
    bp['boxes'][1].set_facecolor('blue')
    plt.ylabel('dN/dS Ratio')
    plt.axhline(y=1, color='gray', linestyle='--', alpha=0.5)
    plt.title('Selection Pressure Comparison')
    
    # Histogram
    plt.subplot(1, 2, 2)
    plt.hist(core_dnds, bins=50, alpha=0.5, label='Core genes', color='blue', density=True)
    plt.hist(operon_dnds, bins=20, alpha=0.7, label='Operon genes', color='red', density=True)
    plt.xlabel('dN/dS Ratio')
    plt.ylabel('Density')
    plt.legend()
    plt.xlim(0, 3)
    plt.axvline(x=1, color='gray', linestyle='--', alpha=0.5)
    plt.title('dN/dS Distribution')
    
    plt.tight_layout()
    output_plot = args.output.replace('.txt', '.png')
    plt.savefig(output_plot, dpi=300, bbox_inches='tight')
    
    print(f"Summary written to: {args.output}")
    print(f"Plot saved to: {output_plot}")

if __name__ == '__main__':
    main()