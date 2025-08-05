#!/usr/bin/env python3
"""
Create gap analysis plots for gene and promoter alignments.
Shows gap frequency across alignment positions.
"""

import os
import pandas as pd
import numpy as np
from Bio import AlignIO
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import argparse

def analyze_gaps_in_alignment(alignment_file):
    """Analyze gap patterns in a multiple sequence alignment."""
    
    if not os.path.exists(alignment_file):
        print(f"Warning: Alignment file {alignment_file} not found")
        return None
    
    try:
        alignment = AlignIO.read(alignment_file, "fasta")
        n_seqs = len(alignment)
        alignment_length = alignment.get_alignment_length()
        
        gap_counts = []
        gap_frequencies = []
        
        for pos in range(alignment_length):
            # Get column
            column = alignment[:, pos]
            
            # Count gaps
            gaps = sum(1 for char in column if char == '-')
            gap_counts.append(gaps)
            gap_frequencies.append(gaps / n_seqs)
        
        return {
            'gap_counts': gap_counts,
            'gap_frequencies': gap_frequencies,
            'n_sequences': n_seqs,
            'alignment_length': alignment_length,
            'total_gaps': sum(gap_counts),
            'avg_gap_frequency': np.mean(gap_frequencies)
        }
    
    except Exception as e:
        print(f"Error analyzing {alignment_file}: {e}")
        return None

def create_gap_analysis_plots(output_dir):
    """Create gap analysis plots for all available alignments."""
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Analyze gene alignments
    gene_dir = "output/msa/dna_alignments"
    gene_files = []
    if os.path.exists(gene_dir):
        gene_files = [f for f in os.listdir(gene_dir) if f.endswith('_aligned.fasta')]
    
    # Analyze promoter alignment (if it exists)
    promoter_file = "output/msa/noncoding_alignments/promoter_aligned.fasta"
    
    all_gap_data = {}
    
    # Process gene alignments
    for gene_file in gene_files:
        gene_name = gene_file.replace('_aligned.fasta', '')
        gene_path = os.path.join(gene_dir, gene_file)
        
        print(f"Analyzing gaps in {gene_name}...")
        gap_data = analyze_gaps_in_alignment(gene_path)
        
        if gap_data:
            gap_data['name'] = gene_name
            gap_data['type'] = 'gene'
            all_gap_data[gene_name] = gap_data
    
    # Process promoter alignment if available
    if os.path.exists(promoter_file):
        print("Analyzing gaps in promoter...")
        gap_data = analyze_gaps_in_alignment(promoter_file)
        if gap_data:
            gap_data['name'] = 'promoter'
            gap_data['type'] = 'promoter'
            all_gap_data['promoter'] = gap_data
    
    if not all_gap_data:
        print("No alignment data found for gap analysis")
        return
    
    # Create individual gap profile plots (matching conservation plot style)
    print("Creating individual gap profile plots...")
    
    for name, data in all_gap_data.items():
        fig, ax = plt.subplots(1, 1, figsize=(12, 4))
        
        positions = range(1, len(data['gap_frequencies']) + 1)
        gap_freqs = data['gap_frequencies']
        
        # Plot with appropriate color
        color = 'darkorange' if data['type'] == 'promoter' else 'blue'
        
        ax.plot(positions, gap_freqs, color=color, linewidth=0.8)
        ax.fill_between(positions, gap_freqs, alpha=0.3, color=color)
        
        ax.set_ylim(0, 1)
        ax.set_xlabel('Position')
        ax.set_ylabel('Gap Frequency')
        ax.set_title(f'{name} - Gap Frequency Profile')
        
        # Add average line
        avg_gap_freq = data['avg_gap_frequency']
        ax.axhline(y=avg_gap_freq, color='red', linestyle='--', 
                   label=f'Average: {avg_gap_freq:.3f}')
        ax.legend()
        
        # Add grid
        ax.grid(True, alpha=0.3)
        
        # Add statistics
        max_gap_freq = np.max(gap_freqs)
        ax.text(0.02, 0.98, 
                f'Sequences: {data["n_sequences"]}\n'
                f'Length: {data["alignment_length"]} bp\n'
                f'Max gap freq: {max_gap_freq:.3f}\n'
                f'Total gaps: {data["total_gaps"]}', 
                transform=ax.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        
        # Save individual plot
        output_file = os.path.join(output_dir, f"{name}_gap_profile.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"  Saved {name} gap profile to {output_file}")
    
    # Create combined gap summary plot
    print("Creating combined gap summary plot...")
    
    # Prepare data for summary plots
    gene_data = [data for data in all_gap_data.values() if data['type'] == 'gene']
    promoter_data = [data for data in all_gap_data.values() if data['type'] == 'promoter']
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # 1. Average gap frequency by sequence type
    if gene_data and promoter_data:
        gene_avg_gaps = [data['avg_gap_frequency'] for data in gene_data]
        promoter_avg_gaps = [data['avg_gap_frequency'] for data in promoter_data]
        
        ax1 = axes[0, 0]
        ax1.boxplot([gene_avg_gaps, promoter_avg_gaps], labels=['Genes', 'Promoter'])
        ax1.set_ylabel('Average Gap Frequency')
        ax1.set_title('Gap Frequency by Sequence Type')
    
    # 2. Gap frequency distribution for genes
    ax2 = axes[0, 1]
    if gene_data:
        gene_names = [data['name'] for data in gene_data]
        gene_gap_freqs = [data['avg_gap_frequency'] for data in gene_data]
        
        bars = ax2.bar(gene_names, gene_gap_freqs, color='skyblue', alpha=0.7)
        ax2.set_ylabel('Average Gap Frequency')
        ax2.set_title('Gap Frequency by Gene')
        ax2.tick_params(axis='x', rotation=45)
        
        # Add promoter as comparison if available
        if promoter_data:
            promoter_gap_freq = promoter_data[0]['avg_gap_frequency']
            ax2.axhline(y=promoter_gap_freq, color='orange', linestyle='--', 
                       label=f'Promoter: {promoter_gap_freq:.3f}')
            ax2.legend()
    
    # 3. Total gaps vs sequence length
    ax3 = axes[1, 0]
    if all_gap_data:
        lengths = [data['alignment_length'] for data in all_gap_data.values()]
        total_gaps = [data['total_gaps'] for data in all_gap_data.values()]
        colors = ['orange' if data['type'] == 'promoter' else 'blue' for data in all_gap_data.values()]
        
        scatter = ax3.scatter(lengths, total_gaps, c=colors, alpha=0.7)
        ax3.set_xlabel('Alignment Length (bp)')
        ax3.set_ylabel('Total Gaps')
        ax3.set_title('Total Gaps vs Alignment Length')
        
        # Add labels
        for i, (name, data) in enumerate(all_gap_data.items()):
            ax3.annotate(name, (lengths[i], total_gaps[i]), 
                        xytext=(5, 5), textcoords='offset points', fontsize=8)
    
    # 4. Gap statistics summary
    ax4 = axes[1, 1]
    ax4.axis('off')
    
    # Create summary text
    summary_text = "GAP ANALYSIS SUMMARY\n" + "="*30 + "\n\n"
    
    if gene_data:
        gene_avg = np.mean([data['avg_gap_frequency'] for data in gene_data])
        summary_text += f"GENES ({len(gene_data)} analyzed):\n"
        summary_text += f"  Avg gap frequency: {gene_avg:.4f}\n"
        summary_text += f"  Range: {min([data['avg_gap_frequency'] for data in gene_data]):.4f} - "
        summary_text += f"{max([data['avg_gap_frequency'] for data in gene_data]):.4f}\n\n"
        
        # Most and least gappy genes
        most_gappy = max(gene_data, key=lambda x: x['avg_gap_frequency'])
        least_gappy = min(gene_data, key=lambda x: x['avg_gap_frequency'])
        summary_text += f"  Most gappy: {most_gappy['name']} ({most_gappy['avg_gap_frequency']:.4f})\n"
        summary_text += f"  Least gappy: {least_gappy['name']} ({least_gappy['avg_gap_frequency']:.4f})\n\n"
    
    if promoter_data:
        promoter_gap_freq = promoter_data[0]['avg_gap_frequency']
        summary_text += f"PROMOTER:\n"
        summary_text += f"  Gap frequency: {promoter_gap_freq:.4f}\n"
        summary_text += f"  Length: {promoter_data[0]['alignment_length']} bp\n"
        summary_text += f"  Total gaps: {promoter_data[0]['total_gaps']}\n"
    
    ax4.text(0.05, 0.95, summary_text, transform=ax4.transAxes, 
             verticalalignment='top', fontsize=10, fontfamily='monospace')
    
    plt.tight_layout()
    
    # Save combined plot
    combined_file = os.path.join(output_dir, "gap_analysis_summary.png")
    plt.savefig(combined_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved combined gap analysis to {combined_file}")
    
    # Print summary
    print(f"\n" + "="*50)
    print("GAP ANALYSIS SUMMARY")
    print("="*50)
    
    for name, data in all_gap_data.items():
        print(f"\n{name.upper()} ({data['type']}):")
        print(f"  Sequences: {data['n_sequences']}")
        print(f"  Length: {data['alignment_length']} bp")
        print(f"  Total gaps: {data['total_gaps']}")
        print(f"  Average gap frequency: {data['avg_gap_frequency']:.4f}")
        print(f"  Max gap frequency: {max(data['gap_frequencies']):.4f}")

def move_promoter_plots_to_diversity_dir():
    """Move promoter plots to the diversity_analysis directory."""
    
    source_dir = "output/promoter_analysis"
    target_dir = "output/diversity_analysis"
    
    if not os.path.exists(source_dir):
        print("No promoter analysis directory found")
        return
    
    print("Moving promoter plots to diversity_analysis directory...")
    
    promoter_files = [
        "promoter_conservation_profile.png",
        "promoter_conservation_summary.png"
    ]
    
    for filename in promoter_files:
        source_file = os.path.join(source_dir, filename)
        target_file = os.path.join(target_dir, filename)
        
        if os.path.exists(source_file):
            import shutil
            shutil.copy2(source_file, target_file)
            print(f"  Copied {filename} to {target_dir}")

def main():
    parser = argparse.ArgumentParser(description="Create gap analysis plots")
    parser.add_argument("--output-dir", default="output/diversity_analysis",
                        help="Output directory for plots")
    
    args = parser.parse_args()
    
    print("Creating gap analysis plots...")
    print("="*50)
    
    # Create gap analysis plots
    create_gap_analysis_plots(args.output_dir)
    
    # Move promoter plots to same directory
    move_promoter_plots_to_diversity_dir()
    
    print(f"\n✅ All plots saved to: {args.output_dir}")
    print("Files created:")
    print("  - Individual gap profiles: *_gap_profile.png")
    print("  - Combined gap summary: gap_analysis_summary.png")
    print("  - Promoter conservation: promoter_conservation_profile.png")

if __name__ == "__main__":
    main()