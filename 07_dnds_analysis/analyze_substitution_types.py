#!/usr/bin/env python3
"""Analyze substitution types in MSAs - synonymous vs non-synonymous"""

import os
import glob
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from collections import defaultdict
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

def analyze_msa_substitutions(msa_file, reference_seq=None):
    """Analyze substitutions in an MSA file"""
    
    try:
        alignment = AlignIO.read(msa_file, 'fasta')
    except:
        print(f"Error reading {msa_file}")
        return None
    
    results = {
        'total_sequences': len(alignment),
        'alignment_length': alignment.get_alignment_length(),
        'synonymous': 0,
        'non_synonymous': 0,
        'indels': 0,
        'total_substitutions': 0
    }
    
    # If no reference provided, use first sequence
    if reference_seq is None:
        reference = str(alignment[0].seq).upper()
    else:
        reference = str(reference_seq).upper()
    
    # Analyze each sequence against reference
    for record in alignment[1:]:  # Skip first if it's the reference
        seq = str(record.seq).upper()
        
        # Compare in codon triplets
        for i in range(0, len(reference) - 2, 3):
            ref_codon = reference[i:i+3]
            seq_codon = seq[i:i+3]
            
            if ref_codon == seq_codon:
                continue
            
            # Check for indels
            if '-' in ref_codon or '-' in seq_codon:
                results['indels'] += 1
                continue
            
            # Classify substitution
            change_type = get_codon_changes(ref_codon, seq_codon)
            if change_type == 'synonymous':
                results['synonymous'] += 1
            elif change_type == 'non-synonymous':
                results['non_synonymous'] += 1
            
            results['total_substitutions'] += 1
    
    # Calculate ratios
    if results['synonymous'] > 0:
        results['dn_ds_ratio'] = results['non_synonymous'] / results['synonymous']
    else:
        results['dn_ds_ratio'] = float('inf') if results['non_synonymous'] > 0 else 0
    
    return results

def main():
    print("=== Analyzing Substitution Types ===\n")
    
    # 1. Analyze operon genes
    print("1. Analyzing operon genes...")
    operon_results = {}
    
    # Read reference sequences
    operon_refs = {}
    for record in SeqIO.parse('operon_genes.fasta', 'fasta'):
        gene_name = record.id.split('|')[0]
        operon_refs[gene_name] = record.seq
    
    # Analyze each operon MSA
    for msa_file in glob.glob('msa_output/*.aln'):
        gene_name = os.path.basename(msa_file).replace('.aln', '')
        
        # Find matching reference
        ref_seq = None
        for ref_name, seq in operon_refs.items():
            if gene_name in ref_name or ref_name in gene_name:
                ref_seq = seq
                break
        
        results = analyze_msa_substitutions(msa_file, ref_seq)
        if results:
            operon_results[gene_name] = results
            print(f"  {gene_name}: dN/dS = {results['dn_ds_ratio']:.3f}")
    
    # 2. Analyze sample of core genes
    print("\n2. Analyzing core genes (sample)...")
    core_results = {}
    
    # Sample some core gene alignments
    core_gene_files = glob.glob('core_gene_sequences_95pct/*.fasta')[:50]  # Sample 50
    
    for fasta_file in core_gene_files:
        gene_name = os.path.basename(fasta_file).replace('.fasta', '')
        
        # For core genes, we'll need to align them first
        # Check if alignment exists
        aligned_file = fasta_file.replace('.fasta', '.aln')
        
        if not os.path.exists(aligned_file):
            # Run quick alignment
            cmd = f"mafft --retree 1 --maxiterate 0 {fasta_file} > {aligned_file} 2>/dev/null"
            os.system(cmd)
        
        if os.path.exists(aligned_file):
            results = analyze_msa_substitutions(aligned_file)
            if results:
                core_results[gene_name] = results
    
    print(f"\nAnalyzed {len(core_results)} core genes")
    
    # 3. Create comparison plots
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Plot 1: dN/dS ratios comparison
    ax1 = axes[0, 0]
    operon_dnds = [r['dn_ds_ratio'] for r in operon_results.values() if r['dn_ds_ratio'] < 5]  # Cap at 5 for visualization
    core_dnds = [r['dn_ds_ratio'] for r in core_results.values() if r['dn_ds_ratio'] < 5]
    
    positions = [1, 2]
    bp = ax1.boxplot([operon_dnds, core_dnds], positions=positions, widths=0.6, patch_artist=True)
    bp['boxes'][0].set_facecolor('red')
    bp['boxes'][1].set_facecolor('blue')
    
    ax1.set_xticklabels(['Operon genes', 'Core genes'])
    ax1.set_ylabel('dN/dS ratio', fontsize=12)
    ax1.set_title('Selection Pressure (dN/dS ratios)', fontsize=14, fontweight='bold')
    ax1.axhline(y=1, color='gray', linestyle='--', alpha=0.5)
    ax1.text(1.5, 1.1, 'Neutral evolution', ha='center', fontsize=10, color='gray')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Substitution types breakdown
    ax2 = axes[0, 1]
    
    # Calculate averages
    operon_syn = np.mean([r['synonymous'] for r in operon_results.values()])
    operon_nonsyn = np.mean([r['non_synonymous'] for r in operon_results.values()])
    core_syn = np.mean([r['synonymous'] for r in core_results.values()])
    core_nonsyn = np.mean([r['non_synonymous'] for r in core_results.values()])
    
    x = np.arange(2)
    width = 0.35
    
    ax2.bar(x - width/2, [operon_syn, core_syn], width, label='Synonymous', color='green', alpha=0.7)
    ax2.bar(x + width/2, [operon_nonsyn, core_nonsyn], width, label='Non-synonymous', color='orange', alpha=0.7)
    
    ax2.set_xticks(x)
    ax2.set_xticklabels(['Operon genes', 'Core genes'])
    ax2.set_ylabel('Average substitutions per gene', fontsize=12)
    ax2.set_title('Types of Substitutions', fontsize=14, fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Plot 3: Individual operon genes
    ax3 = axes[1, 0]
    if operon_results:
        genes = list(operon_results.keys())
        dnds_values = [operon_results[g]['dn_ds_ratio'] if operon_results[g]['dn_ds_ratio'] < 5 else 5 for g in genes]
        
        bars = ax3.bar(range(len(genes)), dnds_values, color='red', alpha=0.7)
        ax3.axhline(y=1, color='gray', linestyle='--', alpha=0.5)
        ax3.set_xticks(range(len(genes)))
        ax3.set_xticklabels([g[:20] + '...' if len(g) > 20 else g for g in genes], rotation=45, ha='right')
        ax3.set_ylabel('dN/dS ratio', fontsize=12)
        ax3.set_title('dN/dS for Individual Operon Genes', fontsize=14, fontweight='bold')
        ax3.grid(True, alpha=0.3, axis='y')
    
    # Plot 4: Summary statistics
    ax4 = axes[1, 1]
    ax4.axis('off')
    
    summary_text = "Summary Statistics\n\n"
    summary_text += "Operon genes:\n"
    if operon_results:
        avg_dnds_operon = np.mean([r['dn_ds_ratio'] for r in operon_results.values() if r['dn_ds_ratio'] < 5])
        summary_text += f"  • Average dN/dS: {avg_dnds_operon:.3f}\n"
        summary_text += f"  • Genes analyzed: {len(operon_results)}\n"
        total_syn = sum(r['synonymous'] for r in operon_results.values())
        total_nonsyn = sum(r['non_synonymous'] for r in operon_results.values())
        summary_text += f"  • Total synonymous: {total_syn}\n"
        summary_text += f"  • Total non-synonymous: {total_nonsyn}\n"
    
    summary_text += "\nCore genes (sample):\n"
    if core_results:
        avg_dnds_core = np.mean([r['dn_ds_ratio'] for r in core_results.values() if r['dn_ds_ratio'] < 5])
        summary_text += f"  • Average dN/dS: {avg_dnds_core:.3f}\n"
        summary_text += f"  • Genes analyzed: {len(core_results)}\n"
        total_syn = sum(r['synonymous'] for r in core_results.values())
        total_nonsyn = sum(r['non_synonymous'] for r in core_results.values())
        summary_text += f"  • Total synonymous: {total_syn}\n"
        summary_text += f"  • Total non-synonymous: {total_nonsyn}\n"
    
    summary_text += "\nInterpretation:\n"
    summary_text += "• dN/dS < 1: Purifying selection\n"
    summary_text += "• dN/dS = 1: Neutral evolution\n"
    summary_text += "• dN/dS > 1: Positive selection"
    
    ax4.text(0.1, 0.9, summary_text, transform=ax4.transAxes, fontsize=12,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig('substitution_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Save detailed results
    all_results = []
    for gene, results in operon_results.items():
        results['gene'] = gene
        results['gene_type'] = 'operon'
        all_results.append(results)
    
    for gene, results in core_results.items():
        results['gene'] = gene
        results['gene_type'] = 'core'
        all_results.append(results)
    
    df = pd.DataFrame(all_results)
    df.to_csv('substitution_analysis_results.csv', index=False)
    
    print("\nResults saved:")
    print("  - substitution_analysis.png")
    print("  - substitution_analysis_results.csv")

if __name__ == '__main__':
    main()