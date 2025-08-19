#!/usr/bin/env python3
"""
Create enhanced conservation plots with Shannon entropy and sequence logos.
Improvements over the SNP-based approach:
1. Shannon entropy-based conservation scores (frequency-weighted)
2. Sequence logos showing nucleotide frequencies
3. Both metrics capture biological conservation more accurately
"""

import os
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from collections import defaultdict
import numpy as np
import pandas as pd
try:
    import logomaker
    LOGOMAKER_AVAILABLE = True
except ImportError:
    LOGOMAKER_AVAILABLE = False
    print("Warning: logomaker not available. Install with: pip install logomaker")

def calculate_shannon_entropy_per_position(alignment_file):
    """Calculate Shannon entropy-based conservation scores per position."""
    sequences = []
    for record in SeqIO.parse(alignment_file, "fasta"):
        sequences.append(str(record.seq).upper())
    
    if not sequences:
        return [], [], None
    
    alignment_length = len(sequences[0])
    conservation_scores = []
    gap_frequencies = []
    position_matrices = []
    
    for pos in range(alignment_length):
        # Get all characters at this position
        column = [seq[pos] for seq in sequences if pos < len(seq)]
        
        # Calculate gap frequency
        gap_count = column.count('-')
        gap_freq = gap_count / len(column) if column else 0
        gap_frequencies.append(gap_freq)
        
        # Count nucleotide frequencies (excluding gaps)
        char_counts = defaultdict(int)
        total_non_gap = 0
        for char in column:
            if char != '-':
                char_counts[char] += 1
                total_non_gap += 1
        
        # Create position weight matrix row for this position
        pwm_row = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        if total_non_gap > 0:
            for nucleotide in ['A', 'T', 'G', 'C']:
                pwm_row[nucleotide] = char_counts[nucleotide] / total_non_gap
        position_matrices.append(pwm_row)
        
        # Calculate Shannon entropy
        if total_non_gap > 0:
            entropy = 0
            for count in char_counts.values():
                if count > 0:
                    freq = count / total_non_gap
                    entropy -= freq * np.log2(freq)
            
            # Normalize entropy to 0-1 scale (max entropy for 4 nucleotides is 2)
            max_entropy = 2.0  # log2(4)
            normalized_entropy = entropy / max_entropy if max_entropy > 0 else 0
            
            # Conservation score = 1 - normalized_entropy
            conservation = 1 - normalized_entropy
        else:
            conservation = 0
        
        conservation_scores.append(conservation)
    
    return conservation_scores, gap_frequencies, position_matrices

def create_sequence_logo(position_matrices, output_file, gene_name, window_size=100):
    """Create sequence logo plots in sliding windows."""
    if not LOGOMAKER_AVAILABLE:
        print(f"Skipping sequence logo for {gene_name} - logomaker not available")
        return False
    
    if not position_matrices:
        return False
    
    # Convert to DataFrame
    df = pd.DataFrame(position_matrices)
    
    # Create sequence logos in windows
    num_positions = len(df)
    num_windows = max(1, (num_positions + window_size - 1) // window_size)
    
    if num_positions <= window_size:
        # Single plot for short sequences
        fig, ax = plt.subplots(figsize=(min(20, num_positions * 0.3), 4))
        logomaker.Logo(df, ax=ax, color_scheme='classic')
        ax.set_title(f'{gene_name} Sequence Logo (n={len(position_matrices)} positions)', 
                    fontsize=14, fontweight='bold')
        ax.set_xlabel('Position in Alignment')
        ax.set_ylabel('Information Content (bits)')
    else:
        # Multiple windows for long sequences
        fig, axes = plt.subplots(num_windows, 1, figsize=(20, 4 * num_windows))
        if num_windows == 1:
            axes = [axes]
        
        for i, ax in enumerate(axes):
            start_pos = i * window_size
            end_pos = min((i + 1) * window_size, num_positions)
            window_df = df.iloc[start_pos:end_pos].copy()
            
            # Reset index to start from 0 for each window
            window_df.index = range(start_pos + 1, end_pos + 1)
            
            logomaker.Logo(window_df, ax=ax, color_scheme='classic')
            ax.set_title(f'{gene_name} Sequence Logo - Positions {start_pos + 1}-{end_pos}', 
                        fontsize=12, fontweight='bold')
            ax.set_xlabel('Position in Alignment')
            ax.set_ylabel('Information (bits)')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    return True

def create_enhanced_conservation_plot(alignment_file, output_file, gene_name, title_suffix: str = ""):
    """Create enhanced conservation plot with Shannon entropy."""
    conservation_scores, gap_frequencies, position_matrices = calculate_shannon_entropy_per_position(alignment_file)
    
    if not conservation_scores:
        print(f"Warning: No sequences found in {alignment_file}")
        return False
    
    # Count sequences used
    num_sequences = 0
    for _ in SeqIO.parse(alignment_file, "fasta"):
        num_sequences += 1

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    
    positions = range(1, len(conservation_scores) + 1)
    
    # Top panel: Conservation scores (Shannon entropy-based)
    ax1.plot(positions, conservation_scores, color='blue', linewidth=1)
    ax1.fill_between(positions, conservation_scores, alpha=0.3, color='blue')
    ax1.set_ylabel('Conservation Score\n(Shannon Entropy)', fontsize=12)
    suffix_text = f" — {title_suffix}" if title_suffix else ""
    ax1.set_title(f'{gene_name} Gene Conservation (n={num_sequences}){suffix_text}', 
                  fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, 1)
    
    # Add conservation level annotations
    ax1.axhline(y=0.9, color='green', linestyle='--', alpha=0.5, label='Highly conserved (>0.9)')
    ax1.axhline(y=0.7, color='orange', linestyle='--', alpha=0.5, label='Moderately conserved (>0.7)')
    ax1.legend(loc='upper right', fontsize=10)
    
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
    
    # Create sequence logo if possible
    logo_file = output_file.replace('_conservation.png', '_sequence_logo.png')
    create_sequence_logo(position_matrices, logo_file, gene_name)
    
    return True

def main():
    parser = argparse.ArgumentParser(description="Create enhanced conservation plots with Shannon entropy and sequence logos")
    parser.add_argument("--msa-dir", default="output/msa/dna_alignments", 
                        help="Directory containing MSA files")
    parser.add_argument("--output-dir", default="output/enhanced_conservation_plots", 
                        help="Output directory for conservation plots")
    parser.add_argument("--title-suffix", default="", 
                        help="Suffix text to append in the plot title")
    parser.add_argument("--logo-window-size", type=int, default=100,
                        help="Window size for sequence logo plots (default: 100)")
    
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
        output_file = os.path.join(args.output_dir, f"{gene_name}_enhanced_conservation.png")
        
        print(f"  Processing {gene_name}...")
        
        success = create_enhanced_conservation_plot(msa_path, output_file, gene_name, args.title_suffix)
        if success:
            successful_plots += 1
            print(f"    ✅ Created: {output_file}")
            logo_file = output_file.replace('_enhanced_conservation.png', '_sequence_logo.png')
            if os.path.exists(logo_file):
                print(f"    ✅ Created: {logo_file}")
        else:
            print(f"    ❌ Failed: {msa_file}")
    
    print(f"\n✅ Successfully created {successful_plots}/{len(msa_files)} enhanced conservation plots")
    print(f"📁 Output directory: {args.output_dir}")
    
    if not LOGOMAKER_AVAILABLE:
        print("\n💡 To enable sequence logos, install logomaker:")
        print("   pip install logomaker")

if __name__ == "__main__":
    main()
