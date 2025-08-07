#!/usr/bin/env python3
"""
Create a promoter conservation plot matching the style of the original conservation_profiles.png
"""

import os
import pandas as pd
import numpy as np
from Bio import AlignIO
import matplotlib.pyplot as plt
import argparse
from collections import defaultdict

def calculate_conservation_per_position(alignment_file):
    """Calculate conservation score for each position in the alignment."""
    alignment = AlignIO.read(alignment_file, "fasta")
    n_seqs = len(alignment)
    alignment_length = alignment.get_alignment_length()
    
    conservation_scores = []
    
    for pos in range(alignment_length):
        # Get column
        column = alignment[:, pos]
        
        # Count occurrences of each character
        char_counts = defaultdict(int)
        for char in column:
            if char != '-':
                char_counts[char] += 1
        
        # Calculate conservation (Shannon entropy based)
        if char_counts:
            total = sum(char_counts.values())
            entropy = 0
            for count in char_counts.values():
                if count > 0:
                    p = count / total
                    entropy -= p * np.log2(p)
            
            # Convert entropy to conservation score (0-1)
            max_entropy = np.log2(4)  # For DNA
            conservation = 1 - (entropy / max_entropy)
        else:
            conservation = 0
        
        conservation_scores.append(conservation)
    
    return conservation_scores

def create_promoter_conservation_plot(alignment_file, output_file):
    """Create promoter conservation plot matching the original style."""
    
    if not os.path.exists(alignment_file):
        print(f"Error: Alignment file {alignment_file} not found")
        return None
    
    # Calculate conservation scores
    print("Calculating conservation scores...")
    scores = calculate_conservation_per_position(alignment_file)
    positions = range(1, len(scores) + 1)
    
    # Create plot matching original style
    fig, ax = plt.subplots(1, 1, figsize=(12, 4))
    
    # Plot line with filled area (matching original style)
    ax.plot(positions, scores, color='darkorange', linewidth=0.8)
    ax.fill_between(positions, scores, alpha=0.3, color='orange')
    
    # Set limits and labels
    ax.set_ylim(0, 1)
    ax.set_xlabel('Position')
    ax.set_ylabel('Conservation')
    ax.set_title('Promoter Conservation Profile')
    
    # Add average line (matching original style)
    avg_conservation = np.mean(scores)
    ax.axhline(y=avg_conservation, color='red', linestyle='--', 
               label=f'Average: {avg_conservation:.3f}')
    ax.legend()
    
    # Add grid for better readability
    ax.grid(True, alpha=0.3)
    
    # Add some statistics as text
    min_conservation = np.min(scores)
    max_conservation = np.max(scores)
    ax.text(0.02, 0.98, f'Min: {min_conservation:.3f}\nMax: {max_conservation:.3f}', 
            transform=ax.transAxes, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    
    # Save plot
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Promoter conservation plot saved to {output_file}")
    print(f"Promoter conservation statistics:")
    print(f"  Average: {avg_conservation:.4f}")
    print(f"  Min: {min_conservation:.4f}")
    print(f"  Max: {max_conservation:.4f}")
    print(f"  Length: {len(scores)} positions")
    
    return output_file

def main():
    parser = argparse.ArgumentParser(description="Create promoter conservation plot")
    parser.add_argument("--alignment", 
                        default="output/msa/noncoding_alignments/promoter_aligned.fasta",
                        help="Path to promoter alignment file")
    parser.add_argument("--output", 
                        default="output/enhanced_diversity_analysis/promoter_conservation_profile.png",
                        help="Output plot file")
    parser.add_argument("--output-dir", 
                        help="Output directory (alternative to --output)")
    
    args = parser.parse_args()
    
    # Handle output directory
    if args.output_dir:
        os.makedirs(args.output_dir, exist_ok=True)
        output_file = os.path.join(args.output_dir, "promoter_conservation_profile.png")
    else:
        output_file = args.output
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Check if alignment file exists
    if not os.path.exists(args.alignment):
        print(f"Error: Alignment file {args.alignment} not found")
        print("\nMake sure you have run the promoter extraction and MSA creation steps:")
        print("  python extract_noncoding_sequences.py")
        print("  python create_msa.py --noncoding-only")
        print("\nOr run the full enhanced pipeline:")
        print("  sbatch run_enhanced_diversity_analysis.sh")
        return 1
    
    # Create the plot
    result = create_promoter_conservation_plot(args.alignment, output_file)
    
    if result:
        print(f"\nSuccess! Promoter conservation plot created.")
        print(f"View the plot: {output_file}")
        return 0
    else:
        print("Failed to create promoter conservation plot")
        return 1

if __name__ == "__main__":
    exit(main())