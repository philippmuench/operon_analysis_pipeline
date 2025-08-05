#!/usr/bin/env python3
"""
Analyze MSA diversity and create visualizations showing:
- Number of variants per position
- Gap frequency
- Conservation scores
- Consensus sequence
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.SeqUtils import GC
from collections import Counter
import argparse
from datetime import datetime

def calculate_position_statistics(alignment):
    """Calculate various statistics for each position in the alignment."""
    
    align_length = alignment.get_alignment_length()
    num_seqs = len(alignment)
    
    # Initialize arrays
    variants_per_pos = []
    gap_frequency = []
    conservation_scores = []
    most_common_base = []
    entropy_scores = []
    
    for pos in range(align_length):
        # Get all bases at this position
        column = alignment[:, pos]
        
        # Count occurrences
        base_counts = Counter(column)
        
        # Number of different variants (excluding gaps)
        non_gap_variants = {base: count for base, count in base_counts.items() if base != '-'}
        num_variants = len(non_gap_variants)
        variants_per_pos.append(num_variants)
        
        # Gap frequency
        gap_count = base_counts.get('-', 0)
        gap_freq = gap_count / num_seqs
        gap_frequency.append(gap_freq)
        
        # Conservation score (frequency of most common non-gap base)
        if non_gap_variants:
            max_count = max(non_gap_variants.values())
            conservation = max_count / (num_seqs - gap_count) if (num_seqs - gap_count) > 0 else 0
            most_common = max(non_gap_variants, key=non_gap_variants.get)
        else:
            conservation = 0
            most_common = '-'
        conservation_scores.append(conservation)
        most_common_base.append(most_common)
        
        # Shannon entropy for diversity
        if non_gap_variants:
            total_non_gap = sum(non_gap_variants.values())
            entropy = 0
            for count in non_gap_variants.values():
                if count > 0:
                    freq = count / total_non_gap
                    entropy -= freq * np.log2(freq)
            entropy_scores.append(entropy)
        else:
            entropy_scores.append(0)
    
    return {
        'position': list(range(align_length)),
        'variants': variants_per_pos,
        'gap_frequency': gap_frequency,
        'conservation': conservation_scores,
        'consensus': most_common_base,
        'entropy': entropy_scores
    }

def create_diversity_plots(stats_df, output_prefix, gene_name):
    """Create comprehensive diversity plots."""
    
    fig, axes = plt.subplots(5, 1, figsize=(16, 12), sharex=True)
    fig.suptitle(f'Sequence Diversity Analysis: {gene_name}', fontsize=16)
    
    # Plot 1: Number of variants per position
    ax1 = axes[0]
    ax1.bar(stats_df['position'], stats_df['variants'], width=1, color='darkblue', alpha=0.7)
    ax1.set_ylabel('Number of Variants')
    ax1.set_title('Sequence Variants per Position (excluding gaps)')
    ax1.grid(axis='y', alpha=0.3)
    
    # Plot 2: Gap frequency
    ax2 = axes[1]
    ax2.bar(stats_df['position'], stats_df['gap_frequency'], width=1, color='red', alpha=0.7)
    ax2.set_ylabel('Gap Frequency')
    ax2.set_title('Gap Frequency per Position')
    ax2.set_ylim(0, 1)
    ax2.grid(axis='y', alpha=0.3)
    
    # Plot 3: Conservation score
    ax3 = axes[2]
    ax3.bar(stats_df['position'], stats_df['conservation'], width=1, color='green', alpha=0.7)
    ax3.set_ylabel('Conservation Score')
    ax3.set_title('Conservation Score per Position (1 = fully conserved)')
    ax3.set_ylim(0, 1)
    ax3.grid(axis='y', alpha=0.3)
    
    # Plot 4: Shannon entropy
    ax4 = axes[3]
    ax4.bar(stats_df['position'], stats_df['entropy'], width=1, color='purple', alpha=0.7)
    ax4.set_ylabel('Shannon Entropy')
    ax4.set_title('Sequence Diversity (Shannon Entropy)')
    ax4.grid(axis='y', alpha=0.3)
    
    # Plot 5: Combined view with rolling average
    ax5 = axes[4]
    window_size = 50
    rolling_variants = stats_df['variants'].rolling(window=window_size, center=True).mean()
    rolling_conservation = stats_df['conservation'].rolling(window=window_size, center=True).mean()
    
    ax5_twin = ax5.twinx()
    
    line1 = ax5.plot(stats_df['position'], rolling_variants, 'b-', label=f'Variants ({window_size}bp window)', linewidth=2)
    line2 = ax5_twin.plot(stats_df['position'], rolling_conservation, 'g-', label=f'Conservation ({window_size}bp window)', linewidth=2)
    
    ax5.set_xlabel('Position in Alignment')
    ax5.set_ylabel('Average Number of Variants', color='b')
    ax5_twin.set_ylabel('Average Conservation Score', color='g')
    ax5.set_title(f'Rolling Average (window size: {window_size}bp)')
    ax5.tick_params(axis='y', labelcolor='b')
    ax5_twin.tick_params(axis='y', labelcolor='g')
    
    # Combine legends
    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax5.legend(lines, labels, loc='upper right')
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_diversity_analysis.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_diversity_analysis.pdf', bbox_inches='tight')
    plt.close()

def create_summary_report(alignment, stats_df, output_prefix, gene_name):
    """Create a summary report of the analysis."""
    
    num_seqs = len(alignment)
    align_length = alignment.get_alignment_length()
    
    # Calculate summary statistics
    avg_variants = stats_df['variants'].mean()
    max_variants = stats_df['variants'].max()
    avg_conservation = stats_df['conservation'].mean()
    avg_gap_freq = stats_df['gap_frequency'].mean()
    
    # Find highly variable regions (top 5%)
    threshold_95 = stats_df['variants'].quantile(0.95)
    variable_positions = stats_df[stats_df['variants'] >= threshold_95]['position'].tolist()
    
    # Find highly conserved regions (conservation > 0.95)
    conserved_positions = stats_df[stats_df['conservation'] > 0.95]['position'].tolist()
    
    report = f"""
# MSA Diversity Analysis Report
## Gene: {gene_name}
## Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

### Alignment Statistics
- Number of sequences: {num_seqs:,}
- Alignment length: {align_length:,} bp
- Average variants per position: {avg_variants:.2f}
- Maximum variants at any position: {max_variants}
- Average conservation score: {avg_conservation:.3f}
- Average gap frequency: {avg_gap_freq:.3f}

### Highly Variable Regions (top 5%)
- Number of positions: {len(variable_positions)}
- Threshold (≥ variants): {threshold_95:.0f}
- Position ranges: {_format_ranges(variable_positions)}

### Highly Conserved Regions (>95% conservation)
- Number of positions: {len(conserved_positions)}
- Position ranges: {_format_ranges(conserved_positions)}

### Position-specific Statistics
- See accompanying CSV file for detailed per-position data
- Visualization saved as PNG and PDF files
"""
    
    with open(f'{output_prefix}_summary_report.txt', 'w') as f:
        f.write(report)
    
    # Save detailed statistics
    stats_df.to_csv(f'{output_prefix}_position_statistics.csv', index=False)
    
    return report

def _format_ranges(positions):
    """Format a list of positions into ranges."""
    if not positions:
        return "None"
    
    ranges = []
    start = positions[0]
    end = positions[0]
    
    for pos in positions[1:]:
        if pos == end + 1:
            end = pos
        else:
            if start == end:
                ranges.append(str(start))
            else:
                ranges.append(f"{start}-{end}")
            start = pos
            end = pos
    
    # Add the last range
    if start == end:
        ranges.append(str(start))
    else:
        ranges.append(f"{start}-{end}")
    
    # Limit to first 10 ranges for readability
    if len(ranges) > 10:
        return ', '.join(ranges[:10]) + f"... and {len(ranges)-10} more"
    return ', '.join(ranges)

def main():
    parser = argparse.ArgumentParser(description='Analyze MSA diversity and create visualizations')
    parser.add_argument('alignment_file', help='Input MSA file (FASTA format)')
    parser.add_argument('-o', '--output-prefix', help='Output file prefix', default='msa_diversity')
    parser.add_argument('-n', '--name', help='Gene name for labels', default='Gene')
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.alignment_file):
        print(f"Error: Input file '{args.alignment_file}' not found")
        sys.exit(1)
    
    print(f"Reading alignment from {args.alignment_file}...")
    try:
        alignment = AlignIO.read(args.alignment_file, 'fasta')
    except Exception as e:
        print(f"Error reading alignment: {e}")
        sys.exit(1)
    
    print(f"Alignment loaded: {len(alignment)} sequences, {alignment.get_alignment_length()} positions")
    
    # Calculate position statistics
    print("Calculating position-specific statistics...")
    stats = calculate_position_statistics(alignment)
    stats_df = pd.DataFrame(stats)
    
    # Create visualizations
    print("Creating diversity plots...")
    create_diversity_plots(stats_df, args.output_prefix, args.name)
    
    # Create summary report
    print("Generating summary report...")
    report = create_summary_report(alignment, stats_df, args.output_prefix, args.name)
    print(report)
    
    print(f"\nAnalysis complete! Output files:")
    print(f"  - {args.output_prefix}_diversity_analysis.png")
    print(f"  - {args.output_prefix}_diversity_analysis.pdf")
    print(f"  - {args.output_prefix}_position_statistics.csv")
    print(f"  - {args.output_prefix}_summary_report.txt")

if __name__ == "__main__":
    main()