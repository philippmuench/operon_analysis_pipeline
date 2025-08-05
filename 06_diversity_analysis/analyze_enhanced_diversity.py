#!/usr/bin/env python3
"""
Enhanced diversity analysis for both coding and non-coding operon sequences.
Calculates conservation metrics, substitution patterns, and generates visualizations.
"""

import os
import pandas as pd
import numpy as np
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import argparse

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

def calculate_pairwise_identity(alignment_file):
    """Calculate average pairwise sequence identity."""
    alignment = AlignIO.read(alignment_file, "fasta")
    n_seqs = len(alignment)
    
    if n_seqs < 2:
        return 0.0
    
    identities = []
    for i in range(n_seqs):
        for j in range(i + 1, n_seqs):
            seq1 = str(alignment[i].seq).replace('-', '')
            seq2 = str(alignment[j].seq).replace('-', '')
            
            if len(seq1) == 0 or len(seq2) == 0:
                continue
            
            # Calculate identity for overlapping region
            min_len = min(len(seq1), len(seq2))
            matches = sum(1 for k in range(min_len) if seq1[k] == seq2[k])
            identity = matches / min_len * 100
            identities.append(identity)
    
    return np.mean(identities) if identities else 0.0

def count_substitutions(alignment_file):
    """Count substitution patterns in the alignment."""
    alignment = AlignIO.read(alignment_file, "fasta")
    
    transitions = ['AG', 'GA', 'CT', 'TC']
    transversions = ['AC', 'CA', 'AT', 'TA', 'GC', 'CG', 'GT', 'TG']
    
    ti_count = 0
    tv_count = 0
    total_subs = 0
    
    # Compare first sequence with all others as reference
    if len(alignment) < 2:
        return 0, 0, 0
    
    ref_seq = str(alignment[0].seq).replace('-', '').upper()
    
    for i in range(1, len(alignment)):
        comp_seq = str(alignment[i].seq).replace('-', '').upper()
        min_len = min(len(ref_seq), len(comp_seq))
        
        for pos in range(min_len):
            ref_base = ref_seq[pos]
            comp_base = comp_seq[pos]
            
            if ref_base != comp_base and ref_base in 'ATGC' and comp_base in 'ATGC':
                sub_type = ref_base + comp_base
                total_subs += 1
                
                if sub_type in transitions:
                    ti_count += 1
                elif sub_type in transversions:
                    tv_count += 1
    
    return ti_count, tv_count, total_subs

def analyze_sequence_diversity(alignment_dir, sequence_type="coding"):
    """Analyze diversity for all alignments in a directory."""
    
    if not os.path.exists(alignment_dir):
        print(f"Warning: Alignment directory {alignment_dir} not found")
        return pd.DataFrame()
    
    alignment_files = [f for f in os.listdir(alignment_dir) if f.endswith('_aligned.fasta')]
    
    if not alignment_files:
        print(f"Warning: No alignment files found in {alignment_dir}")
        return pd.DataFrame()
    
    results = []
    
    for align_file in alignment_files:
        gene_name = align_file.replace('_aligned.fasta', '')
        align_path = os.path.join(alignment_dir, align_file)
        
        print(f"Analyzing {sequence_type} diversity for {gene_name}...")
        
        try:
            # Load alignment
            alignment = AlignIO.read(align_path, "fasta")
            n_sequences = len(alignment)
            alignment_length = alignment.get_alignment_length()
            
            # Calculate conservation
            conservation_scores = calculate_conservation_per_position(align_path)
            avg_conservation = np.mean(conservation_scores)
            min_conservation = np.min(conservation_scores)
            
            # Calculate pairwise identity
            pairwise_identity = calculate_pairwise_identity(align_path)
            
            # Count substitutions
            transitions, transversions, total = count_substitutions(align_path)
            
            # Calculate Ti/Tv ratio
            ti_tv_ratio = transitions / transversions if transversions > 0 else np.inf
            
            result = {
                'sequence_name': gene_name,
                'sequence_type': sequence_type,
                'n_sequences': n_sequences,
                'alignment_length': alignment_length,
                'avg_conservation': avg_conservation,
                'min_conservation': min_conservation,
                'pairwise_identity': pairwise_identity,
                'transitions': transitions,
                'transversions': transversions,
                'total_substitutions': total,
                'ti_tv_ratio': ti_tv_ratio,
                'conservation_scores': conservation_scores  # Store for plotting
            }
            
            results.append(result)
            
        except Exception as e:
            print(f"Error analyzing {gene_name}: {e}")
    
    return pd.DataFrame(results)

def create_conservation_plots(results_df, output_dir):
    """Create conservation profile plots."""
    
    # Separate coding and non-coding
    coding_results = results_df[results_df['sequence_type'] == 'coding']
    noncoding_results = results_df[results_df['sequence_type'] == 'noncoding']
    
    # Plot coding sequences
    if not coding_results.empty:
        fig, axes = plt.subplots(len(coding_results), 1, figsize=(12, 3*len(coding_results)))
        if len(coding_results) == 1:
            axes = [axes]
        
        for i, (_, row) in enumerate(coding_results.iterrows()):
            axes[i].plot(row['conservation_scores'], linewidth=1)
            axes[i].set_title(f"{row['sequence_name']} - Conservation Profile")
            axes[i].set_ylabel("Conservation Score")
            axes[i].set_ylim(0, 1)
            axes[i].grid(True, alpha=0.3)
        
        axes[-1].set_xlabel("Alignment Position")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "coding_conservation_profiles.png"), dpi=300, bbox_inches='tight')
        plt.close()
    
    # Plot non-coding sequences
    if not noncoding_results.empty:
        fig, axes = plt.subplots(len(noncoding_results), 1, figsize=(12, 3*len(noncoding_results)))
        if len(noncoding_results) == 1:
            axes = [axes]
        
        for i, (_, row) in enumerate(noncoding_results.iterrows()):
            axes[i].plot(row['conservation_scores'], linewidth=1, color='orange')
            axes[i].set_title(f"{row['sequence_name']} - Conservation Profile")
            axes[i].set_ylabel("Conservation Score")
            axes[i].set_ylim(0, 1)
            axes[i].grid(True, alpha=0.3)
        
        axes[-1].set_xlabel("Alignment Position")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "noncoding_conservation_profiles.png"), dpi=300, bbox_inches='tight')
        plt.close()

def create_promoter_conservation_plot(results_df, output_dir):
    """Create promoter conservation plot matching the original style."""
    
    # Get promoter data
    promoter_results = results_df[
        (results_df['sequence_type'] == 'noncoding') & 
        (results_df['sequence_name'] == 'promoter')
    ]
    
    if promoter_results.empty:
        print("No promoter data found for conservation plot")
        return
    
    # Create plot matching original style
    fig, ax = plt.subplots(1, 1, figsize=(12, 4))
    
    promoter_row = promoter_results.iloc[0]
    scores = promoter_row['conservation_scores']
    positions = range(1, len(scores) + 1)
    
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
    
    plt.tight_layout()
    
    # Save with same naming convention as original
    output_file = os.path.join(output_dir, "promoter_conservation_profile.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Promoter conservation plot saved to {output_file}")
    
    return output_file

def create_summary_plot(results_df, output_dir):
    """Create summary comparison plots."""
    
    if results_df.empty:
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Conservation comparison
    sns.boxplot(data=results_df, x='sequence_type', y='avg_conservation', ax=axes[0,0])
    axes[0,0].set_title('Average Conservation by Sequence Type')
    axes[0,0].set_ylabel('Conservation Score')
    
    # Pairwise identity comparison
    sns.boxplot(data=results_df, x='sequence_type', y='pairwise_identity', ax=axes[0,1])
    axes[0,1].set_title('Pairwise Identity by Sequence Type')
    axes[0,1].set_ylabel('Pairwise Identity (%)')
    
    # Ti/Tv ratio comparison
    # Filter out infinite values for plotting
    plot_data = results_df[results_df['ti_tv_ratio'] != np.inf].copy()
    if not plot_data.empty:
        sns.boxplot(data=plot_data, x='sequence_type', y='ti_tv_ratio', ax=axes[1,0])
        axes[1,0].set_title('Ti/Tv Ratio by Sequence Type')
        axes[1,0].set_ylabel('Ti/Tv Ratio')
    
    # Conservation vs Identity scatter
    sns.scatterplot(data=results_df, x='avg_conservation', y='pairwise_identity', 
                   hue='sequence_type', ax=axes[1,1])
    axes[1,1].set_title('Conservation vs Pairwise Identity')
    axes[1,1].set_xlabel('Average Conservation')
    axes[1,1].set_ylabel('Pairwise Identity (%)')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "enhanced_diversity_summary.png"), dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Enhanced diversity analysis for coding and non-coding sequences")
    parser.add_argument("--coding-alignments", default="output/msa/dna_alignments", 
                        help="Directory containing coding sequence alignments")
    parser.add_argument("--noncoding-alignments", default="output/msa/noncoding_alignments", 
                        help="Directory containing non-coding sequence alignments")
    parser.add_argument("--output-dir", default="output/enhanced_diversity_analysis", 
                        help="Output directory for results")
    parser.add_argument("--coding-only", action="store_true", 
                        help="Only analyze coding sequences")
    parser.add_argument("--noncoding-only", action="store_true", 
                        help="Only analyze non-coding sequences")
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    all_results = []
    
    # Analyze coding sequences
    if not args.noncoding_only:
        print("Analyzing coding sequence diversity...")
        coding_results = analyze_sequence_diversity(args.coding_alignments, "coding")
        if not coding_results.empty:
            all_results.append(coding_results)
    
    # Analyze non-coding sequences
    if not args.coding_only:
        print("Analyzing non-coding sequence diversity...")
        noncoding_results = analyze_sequence_diversity(args.noncoding_alignments, "noncoding")
        if not noncoding_results.empty:
            all_results.append(noncoding_results)
    
    if not all_results:
        print("No results to analyze")
        return
    
    # Combine results
    combined_results = pd.concat(all_results, ignore_index=True)
    
    # Save detailed results (excluding conservation_scores for CSV)
    csv_results = combined_results.drop('conservation_scores', axis=1)
    csv_file = os.path.join(args.output_dir, "enhanced_diversity_results.csv")
    csv_results.to_csv(csv_file, index=False)
    print(f"Detailed results saved to {csv_file}")
    
    # Create visualizations
    print("Creating conservation plots...")
    create_conservation_plots(combined_results, args.output_dir)
    
    # Create promoter-specific plot matching original style
    print("Creating promoter conservation plot...")
    create_promoter_conservation_plot(combined_results, args.output_dir)
    
    print("Creating summary plots...")
    create_summary_plot(combined_results, args.output_dir)
    
    # Print summary statistics
    print("\n" + "="*60)
    print("ENHANCED DIVERSITY ANALYSIS SUMMARY")
    print("="*60)
    
    for seq_type in ['coding', 'noncoding']:
        type_results = combined_results[combined_results['sequence_type'] == seq_type]
        if not type_results.empty:
            print(f"\n{seq_type.upper()} SEQUENCES:")
            print(f"Number of sequences analyzed: {len(type_results)}")
            print(f"Average conservation: {type_results['avg_conservation'].mean():.4f}")
            print(f"Average pairwise identity: {type_results['pairwise_identity'].mean():.2f}%")
            print(f"Average Ti/Tv ratio: {type_results[type_results['ti_tv_ratio'] != np.inf]['ti_tv_ratio'].mean():.2f}")
            
            # Most and least conserved
            most_conserved = type_results.loc[type_results['avg_conservation'].idxmax()]
            least_conserved = type_results.loc[type_results['avg_conservation'].idxmin()]
            print(f"Most conserved: {most_conserved['sequence_name']} ({most_conserved['avg_conservation']:.4f})")
            print(f"Least conserved: {least_conserved['sequence_name']} ({least_conserved['avg_conservation']:.4f})")
    
    print(f"\nResults saved to: {args.output_dir}")

if __name__ == "__main__":
    main()