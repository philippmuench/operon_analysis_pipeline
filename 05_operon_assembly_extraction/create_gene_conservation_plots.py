#!/usr/bin/env python3
"""
Create SNP-based conservation plots for all operon genes.
Similar to the promoter plot but for coding sequences.
"""

import os
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from collections import defaultdict
import numpy as np

def calculate_conservation_per_position(alignment_file):
    """Calculate SNP-based conservation and gap frequencies per position."""
    sequences = []
    for record in SeqIO.parse(alignment_file, "fasta"):
        sequences.append(str(record.seq).upper())
    
    if not sequences:
        return [], []
    
    alignment_length = len(sequences[0])
    conservation_scores = []
    gap_frequencies = []
    
    for pos in range(alignment_length):
        # Get all characters at this position
        column = [seq[pos] for seq in sequences if pos < len(seq)]
        
        # Calculate gap frequency
        gap_count = column.count('-')
        gap_freq = gap_count / len(column) if column else 0
        gap_frequencies.append(gap_freq)
        
        # Calculate conservation (SNP-based, excluding gaps)
        char_counts = defaultdict(int)
        for char in column:
            if char != '-':
                char_counts[char] += 1
        
        # Calculate SNP count (excluding gaps)
        if char_counts:
            num_variants = len(char_counts)
            snp_count = max(0, num_variants - 1)
            
            # Maximum possible SNPs is 3 (A, T, G, C variants)
            max_possible_snps = 3
            conservation = 1 - (snp_count / max_possible_snps)
        else:
            conservation = 0
        
        conservation_scores.append(conservation)
    
    return conservation_scores, gap_frequencies

def create_gene_conservation_plot(alignment_file, output_file, gene_name):
    """Create a two-panel conservation plot for a gene."""
    conservation_scores, gap_frequencies = calculate_conservation_per_position(alignment_file)
    
    if not conservation_scores:
        print(f"Warning: No sequences found in {alignment_file}")
        return False
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    
    positions = range(1, len(conservation_scores) + 1)
    
    # Top panel: Conservation scores
    ax1.plot(positions, conservation_scores, color='blue', linewidth=1)
    ax1.fill_between(positions, conservation_scores, alpha=0.3, color='blue')
    ax1.set_ylabel('Conservation Score\n(SNP-based)', fontsize=12)
    ax1.set_title(f'{gene_name} Gene Conservation Analysis', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, 1)
    
    # Bottom panel: Gap frequencies
    ax2.plot(positions, gap_frequencies, color='red', linewidth=1)
    ax2.fill_between(positions, gap_frequencies, alpha=0.3, color='red')
    ax2.set_ylabel('Gap Frequency', fontsize=12)
    ax2.set_xlabel('Position in Alignment', fontsize=12)
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0, 1)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    return True

def main():
    parser = argparse.ArgumentParser(description="Create conservation plots for operon genes")
    parser.add_argument("--msa-dir", default="output/msa/dna_alignments", 
                        help="Directory containing MSA files")
    parser.add_argument("--output-dir", default="output/plots/gene_conservation", 
                        help="Output directory for conservation plots")
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Process all MSA files
    msa_files = [f for f in os.listdir(args.msa_dir) if f.endswith('_aligned.fasta')]
    
    if not msa_files:
        print(f"No MSA files found in {args.msa_dir}")
        return
    
    print(f"Found {len(msa_files)} MSA files to process:")
    
    successful_plots = 0
    for msa_file in sorted(msa_files):
        gene_name = msa_file.replace('_aligned.fasta', '')
        msa_path = os.path.join(args.msa_dir, msa_file)
        output_file = os.path.join(args.output_dir, f"{gene_name}_conservation.png")
        
        print(f"  Processing {gene_name}...")
        
        success = create_gene_conservation_plot(msa_path, output_file, gene_name)
        if success:
            successful_plots += 1
            print(f"    ✅ Created: {output_file}")
        else:
            print(f"    ❌ Failed: {msa_file}")
    
    print(f"\n✅ Successfully created {successful_plots}/{len(msa_files)} conservation plots")
    print(f"📁 Output directory: {args.output_dir}")

if __name__ == "__main__":
    main()
