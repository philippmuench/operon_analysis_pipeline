#!/usr/bin/env python3
"""Analyze substitution types in operon gene MSAs"""

import os
import glob
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def get_codon_changes(codon1, codon2):
    """Classify codon change as synonymous or non-synonymous"""
    if len(codon1) != 3 or len(codon2) != 3:
        return None
    
    # Skip if contains gaps or ambiguous nucleotides
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

def analyze_msa_pairwise(msa_file, max_pairs=100):
    """Analyze substitutions in an MSA file using pairwise comparisons"""
    
    try:
        alignment = AlignIO.read(msa_file, 'fasta')
    except:
        print(f"Error reading {msa_file}")
        return None
    
    n_seqs = len(alignment)
    results = {
        'gene': os.path.basename(msa_file).replace('_variants_aligned.fasta', ''),
        'n_sequences': n_seqs,
        'alignment_length': alignment.get_alignment_length(),
        'synonymous': 0,
        'non_synonymous': 0,
        'total_substitutions': 0,
        'pairs_analyzed': 0
    }
    
    # Sample sequence pairs if too many
    if n_seqs > 20:
        indices = np.random.choice(n_seqs, 20, replace=False)
    else:
        indices = range(n_seqs)
    
    # Analyze pairwise
    for i in range(len(indices)-1):
        for j in range(i+1, min(i+5, len(indices))):  # Limit comparisons
            seq1 = str(alignment[int(indices[i])].seq).upper()
            seq2 = str(alignment[int(indices[j])].seq).upper()
            
            # Compare in codon triplets
            for pos in range(0, len(seq1) - 2, 3):
                codon1 = seq1[pos:pos+3]
                codon2 = seq2[pos:pos+3]
                
                if codon1 == codon2:
                    continue
                
                change_type = get_codon_changes(codon1, codon2)
                if change_type == 'synonymous':
                    results['synonymous'] += 1
                elif change_type == 'non-synonymous':
                    results['non_synonymous'] += 1
                    
                if change_type:
                    results['total_substitutions'] += 1
            
            results['pairs_analyzed'] += 1
    
    # Calculate ratio
    if results['synonymous'] > 0:
        results['dn_ds_ratio'] = results['non_synonymous'] / results['synonymous']
    else:
        results['dn_ds_ratio'] = float('inf') if results['non_synonymous'] > 0 else 0
    
    # Normalize by pairs analyzed
    if results['pairs_analyzed'] > 0:
        results['syn_per_pair'] = results['synonymous'] / results['pairs_analyzed']
        results['nonsyn_per_pair'] = results['non_synonymous'] / results['pairs_analyzed']
    
    return results

def main():
    print("=== Analyzing Operon Gene Substitutions ===\n")
    
    # Get operon MSA files
    msa_files = glob.glob('msa_output/*_aligned.fasta')
    print(f"Found {len(msa_files)} operon MSA files")
    
    results = []
    for msa_file in msa_files:
        print(f"Analyzing {os.path.basename(msa_file)}...")
        result = analyze_msa_pairwise(msa_file)
        if result:
            results.append(result)
            print(f"  Synonymous: {result['synonymous']}, Non-synonymous: {result['non_synonymous']}")
            print(f"  dN/dS ratio: {result['dn_ds_ratio']:.3f}")
    
    # Create visualization
    if results:
        df = pd.DataFrame(results)
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Plot 1: dN/dS ratios
        genes = df['gene'].str.replace('fructoselysine_glucoselysine_PTS_system_', 'PTS_')
        genes = genes.str.replace('sigma-54_dependent_transcriptional_regulator,_gfr_operon_transcriptional_activator', 'Sigma-54')
        
        dnds_values = df['dn_ds_ratio'].replace([np.inf], 3.0)  # Cap infinity at 3
        
        bars = ax1.bar(range(len(genes)), dnds_values, color='red', alpha=0.7)
        ax1.axhline(y=1, color='gray', linestyle='--', alpha=0.5, label='Neutral (dN/dS=1)')
        ax1.set_xticks(range(len(genes)))
        ax1.set_xticklabels(genes, rotation=45, ha='right')
        ax1.set_ylabel('dN/dS ratio', fontsize=12)
        ax1.set_title('Selection Pressure on Operon Genes', fontsize=14, fontweight='bold')
        ax1.legend()
        ax1.grid(True, alpha=0.3, axis='y')
        
        # Add text for interpretation
        ax1.text(0.02, 0.95, 'dN/dS < 1: Purifying selection', transform=ax1.transAxes, 
                fontsize=10, bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5))
        
        # Plot 2: Substitution counts
        x = np.arange(len(genes))
        width = 0.35
        
        ax2.bar(x - width/2, df['syn_per_pair'], width, label='Synonymous', color='green', alpha=0.7)
        ax2.bar(x + width/2, df['nonsyn_per_pair'], width, label='Non-synonymous', color='orange', alpha=0.7)
        
        ax2.set_xticks(x)
        ax2.set_xticklabels(genes, rotation=45, ha='right')
        ax2.set_ylabel('Substitutions per sequence pair', fontsize=12)
        ax2.set_title('Types of Substitutions in Operon Genes', fontsize=14, fontweight='bold')
        ax2.legend()
        ax2.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        plt.savefig('operon_substitution_analysis.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # Save results
        df.to_csv('operon_substitution_results.csv', index=False)
        
        # Print summary
        print("\n=== Summary ===")
        print(f"Average dN/dS ratio: {df['dn_ds_ratio'][df['dn_ds_ratio'] < 10].mean():.3f}")
        print(f"Total synonymous substitutions: {df['synonymous'].sum()}")
        print(f"Total non-synonymous substitutions: {df['non_synonymous'].sum()}")
        print(f"Overall ratio: {df['non_synonymous'].sum() / df['synonymous'].sum():.3f}")
        
        print("\nFiles saved:")
        print("  - operon_substitution_analysis.png")
        print("  - operon_substitution_results.csv")

if __name__ == '__main__':
    main()