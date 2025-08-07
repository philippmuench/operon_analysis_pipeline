#!/usr/bin/env python3
"""
Create gap-aware promoter conservation plot that properly handles alignment gaps.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from Bio import AlignIO
from collections import Counter
import argparse

def calculate_gap_aware_conservation(alignment):
    """Calculate conservation scores that properly account for gaps."""
    
    length = alignment.get_alignment_length()
    conservation_scores = []
    gap_frequencies = []
    
    for i in range(length):
        # Get all characters at this position
        column = [seq[i].upper() for seq in alignment]
        
        # Count frequencies including gaps
        counts = Counter(column)
        total_seqs = len(column)
        
        # Calculate gap frequency
        gap_freq = counts.get('-', 0) / total_seqs
        gap_frequencies.append(gap_freq)
        
        # Calculate conservation for non-gap positions
        non_gap_chars = [char for char in column if char != '-']
        
        if len(non_gap_chars) == 0:
            # All gaps - no conservation
            conservation_scores.append(0.0)
        elif len(non_gap_chars) == 1:
            # Only one sequence - perfect conservation
            conservation_scores.append(1.0)
        else:
            # Shannon entropy for non-gap characters
            non_gap_counts = Counter(non_gap_chars)
            total_non_gap = len(non_gap_chars)
            
            entropy = 0.0
            for count in non_gap_counts.values():
                if count > 0:
                    p = count / total_non_gap
                    entropy -= p * np.log2(p)
            
            # Normalize entropy (0 = perfect conservation, higher = more variable)
            max_entropy = np.log2(min(4, total_non_gap))  # Max 4 nucleotides
            if max_entropy > 0:
                normalized_entropy = entropy / max_entropy
                conservation_score = 1.0 - normalized_entropy
            else:
                conservation_score = 1.0
                
            conservation_scores.append(max(0.0, conservation_score))
    
    return conservation_scores, gap_frequencies

def create_gap_aware_plot(msa_file, output_file, pribnow_start=43, pribnow_end=50):
    """Create conservation plot with proper gap handling."""
    
    # Read alignment
    print(f"Reading alignment from {msa_file}...")
    alignment = AlignIO.read(msa_file, "fasta")
    
    print(f"Alignment: {len(alignment)} sequences, {alignment.get_alignment_length()} positions")
    
    # Calculate gap-aware conservation
    conservation_scores, gap_frequencies = calculate_gap_aware_conservation(alignment)
    
    # Create position array
    positions = np.arange(1, len(conservation_scores) + 1)
    
    # Create plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10), height_ratios=[3, 1])
    
    # Main conservation plot
    ax1.fill_between(positions, conservation_scores, alpha=0.3, color='orange', label='Conservation Score')
    ax1.plot(positions, conservation_scores, color='darkorange', linewidth=2)
    
    # Highlight Pribnow box
    if pribnow_start <= len(conservation_scores) and pribnow_end <= len(conservation_scores):
        pribnow_positions = positions[pribnow_start-1:pribnow_end]
        pribnow_scores = conservation_scores[pribnow_start-1:pribnow_end]
        ax1.fill_between(pribnow_positions, pribnow_scores, alpha=0.8, color='red', 
                        label=f'Pribnow Box ({pribnow_start}-{pribnow_end})')
        ax1.axvspan(pribnow_start, pribnow_end, alpha=0.2, color='red')
    
    # Add average line
    avg_conservation = np.mean(conservation_scores)
    ax1.axhline(y=avg_conservation, color='red', linestyle='--', alpha=0.7, 
               label=f'Average Conservation ({avg_conservation:.3f})')
    
    ax1.set_ylabel('Conservation Score', fontsize=12)
    ax1.set_title('Gap-Aware Promoter Conservation Analysis\n(Properly accounting for alignment gaps)', 
                 fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, 1)
    
    # Gap frequency subplot
    ax2.fill_between(positions, gap_frequencies, alpha=0.5, color='lightblue', label='Gap Frequency')
    ax2.plot(positions, gap_frequencies, color='blue', linewidth=1)
    ax2.set_xlabel('Position in Alignment', fontsize=12)
    ax2.set_ylabel('Gap Frequency', fontsize=10)
    ax2.set_title('Gap Frequency per Position', fontsize=12)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0, 1)
    
    # Highlight Pribnow box in gap plot too
    if pribnow_start <= len(gap_frequencies) and pribnow_end <= len(gap_frequencies):
        ax2.axvspan(pribnow_start, pribnow_end, alpha=0.2, color='red')
    
    plt.tight_layout()
    
    # Save plot
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Print statistics
    print(f"\nGap-Aware Conservation Analysis:")
    print(f"==============================")
    print(f"Average conservation: {avg_conservation:.3f}")
    print(f"Min conservation: {np.min(conservation_scores):.3f}")
    print(f"Max conservation: {np.max(conservation_scores):.3f}")
    print(f"Average gap frequency: {np.mean(gap_frequencies):.3f}")
    print(f"Positions with >50% gaps: {np.sum(np.array(gap_frequencies) > 0.5)}")
    print(f"Positions with >90% gaps: {np.sum(np.array(gap_frequencies) > 0.9)}")
    
    if pribnow_start <= len(conservation_scores) and pribnow_end <= len(conservation_scores):
        pribnow_conservation = np.mean(conservation_scores[pribnow_start-1:pribnow_end])
        pribnow_gaps = np.mean(gap_frequencies[pribnow_start-1:pribnow_end])
        print(f"\nPribnow Box (positions {pribnow_start}-{pribnow_end}):")
        print(f"  Average conservation: {pribnow_conservation:.3f}")
        print(f"  Average gap frequency: {pribnow_gaps:.3f}")
    
    print(f"\nPlot saved to: {output_file}")
    
    return conservation_scores, gap_frequencies

if __name__ == "__main__":
    create_gap_aware_plot("output/msa/noncoding_alignments/promoter_aligned.fasta", 
                         "output/plots/gap_aware_promoter_conservation.png")
