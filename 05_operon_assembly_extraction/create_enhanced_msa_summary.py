#!/usr/bin/env python3
"""
Create enhanced MSA summary with conservation and identity metrics.
This exports the same metrics as step 04 core gene analysis for consistency.
"""

import os
import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO, AlignIO
from collections import defaultdict
import glob

def calculate_shannon_entropy_conservation(alignment_file):
    """Calculate Shannon entropy-based conservation scores per position."""
    try:
        alignment = AlignIO.read(alignment_file, "fasta")
        n_seqs = len(alignment)
        alignment_length = alignment.get_alignment_length()
        
        conservation_scores = []
        
        for pos in range(alignment_length):
            column = [seq[pos] for seq in alignment]
            # Remove gaps for conservation calculation
            non_gap_chars = [char for char in column if char != '-']
            
            if len(non_gap_chars) == 0:
                conservation_scores.append(0.0)
                continue
                
            # Calculate Shannon entropy
            char_counts = {}
            for char in non_gap_chars:
                char_counts[char] = char_counts.get(char, 0) + 1
            
            total_chars = len(non_gap_chars)
            entropy = 0
            for count in char_counts.values():
                p = count / total_chars
                if p > 0:
                    entropy -= p * np.log2(p)
            
            # Convert entropy to conservation score (0-1)
            max_entropy = np.log2(min(4, len(char_counts)))  # 4 for DNA
            if max_entropy > 0:
                conservation = 1 - (entropy / max_entropy)
            else:
                conservation = 1.0
                
            conservation_scores.append(conservation)
        
        return conservation_scores
    except Exception as e:
        print(f"Error calculating conservation for {alignment_file}: {e}")
        return []

def calculate_pairwise_identity(alignment_file):
    """Calculate average pairwise sequence identity."""
    try:
        alignment = AlignIO.read(alignment_file, "fasta")
        n_seqs = len(alignment)
        
        if n_seqs < 2:
            return 100.0, 0.0, 100.0, 100.0
        
        identities = []
        
        for i in range(min(n_seqs, 100)):  # Sample max 100 sequences for efficiency
            for j in range(i + 1, min(n_seqs, 100)):
                seq1 = str(alignment[i].seq)
                seq2 = str(alignment[j].seq)
                
                # Calculate identity (excluding gaps)
                matches = 0
                valid_positions = 0
                
                for k in range(len(seq1)):
                    if seq1[k] != '-' and seq2[k] != '-':
                        valid_positions += 1
                        if seq1[k] == seq2[k]:
                            matches += 1
                
                if valid_positions > 0:
                    identity = (matches / valid_positions) * 100
                    identities.append(identity)
        
        if identities:
            return np.mean(identities), np.std(identities), np.min(identities), np.max(identities)
        else:
            return 100.0, 0.0, 100.0, 100.0
            
    except Exception as e:
        print(f"Error calculating identity for {alignment_file}: {e}")
        return 100.0, 0.0, 100.0, 100.0

def calculate_gap_statistics(alignment_file):
    """Calculate gap statistics."""
    try:
        alignment = AlignIO.read(alignment_file, "fasta")
        n_seqs = len(alignment)
        alignment_length = alignment.get_alignment_length()
        
        total_gaps = 0
        gap_positions = 0
        
        for pos in range(alignment_length):
            column = [seq[pos] for seq in alignment]
            gaps_in_column = column.count('-')
            
            if gaps_in_column > 0:
                gap_positions += 1
                total_gaps += gaps_in_column
        
        gap_percentage = (total_gaps / (n_seqs * alignment_length)) * 100
        mean_gaps_per_position = total_gaps / alignment_length
        
        return gap_percentage, mean_gaps_per_position
        
    except Exception as e:
        print(f"Error calculating gaps for {alignment_file}: {e}")
        return 0.0, 0.0

def analyze_single_msa(msa_file):
    """Analyze a single MSA file and return comprehensive metrics."""
    try:
        # Basic info
        alignment = AlignIO.read(msa_file, "fasta")
        n_sequences = len(alignment)
        alignment_length = alignment.get_alignment_length()
        
        # Gene name from filename
        gene_name = os.path.basename(msa_file).replace('_aligned.fasta', '')
        
        # Calculate conservation
        conservation_scores = calculate_shannon_entropy_conservation(msa_file)
        
        # Calculate identity
        mean_identity, std_identity, min_identity, max_identity = calculate_pairwise_identity(msa_file)
        
        # Calculate gaps
        gap_percentage, mean_gaps_per_position = calculate_gap_statistics(msa_file)
        
        # Conservation statistics
        if conservation_scores:
            mean_conservation = np.mean(conservation_scores)
            std_conservation = np.std(conservation_scores)
            min_conservation = np.min(conservation_scores)
            max_conservation = np.max(conservation_scores)
            
            # Find most/least conserved windows (10bp)
            window_size = 10
            window_conservations = []
            for i in range(len(conservation_scores) - window_size + 1):
                window_cons = np.mean(conservation_scores[i:i+window_size])
                window_conservations.append(window_cons)
            
            if window_conservations:
                most_conserved_window = np.max(window_conservations)
                least_conserved_window = np.min(window_conservations)
            else:
                most_conserved_window = mean_conservation
                least_conserved_window = mean_conservation
        else:
            mean_conservation = std_conservation = 0.0
            min_conservation = max_conservation = 0.0
            most_conserved_window = least_conserved_window = 0.0
        
        # Conservation category
        if mean_conservation >= 0.95:
            conservation_category = "Highly conserved"
        elif mean_conservation >= 0.8:
            conservation_category = "Moderately conserved"
        elif mean_conservation >= 0.6:
            conservation_category = "Variable"
        else:
            conservation_category = "Highly variable"
        
        return {
            'gene': gene_name,
            'n_sequences': n_sequences,
            'alignment_length': alignment_length,
            'mean_conservation': mean_conservation,
            'std_conservation': std_conservation,
            'min_conservation': min_conservation,
            'max_conservation': max_conservation,
            'mean_pairwise_identity': mean_identity / 100,  # Convert to 0-1 scale like step 04
            'std_pairwise_identity': std_identity / 100,
            'min_pairwise_identity': min_identity / 100,
            'max_pairwise_identity': max_identity / 100,
            'gap_percentage': gap_percentage,
            'mean_gaps_per_position': mean_gaps_per_position,
            'most_conserved_window': most_conserved_window,
            'least_conserved_window': least_conserved_window,
            'conservation_category': conservation_category,
            'alignment_file': msa_file
        }
        
    except Exception as e:
        print(f"Error analyzing {msa_file}: {e}")
        return None

def analyze_strategy_msas(msa_dir, strategy_name):
    """Analyze all MSAs for a specific strategy."""
    dna_alignments_dir = os.path.join(msa_dir, 'dna_alignments')
    
    if not os.path.exists(dna_alignments_dir):
        print(f"Warning: DNA alignments directory not found: {dna_alignments_dir}")
        return pd.DataFrame()
    
    msa_files = glob.glob(os.path.join(dna_alignments_dir, '*_aligned.fasta'))
    
    if not msa_files:
        print(f"Warning: No MSA files found in {dna_alignments_dir}")
        return pd.DataFrame()
    
    print(f"Analyzing {len(msa_files)} MSA files for {strategy_name}...")
    
    results = []
    for msa_file in msa_files:
        print(f"  Processing {os.path.basename(msa_file)}...")
        result = analyze_single_msa(msa_file)
        if result:
            results.append(result)
    
    return pd.DataFrame(results)

def main():
    parser = argparse.ArgumentParser(description="Create enhanced MSA summary with conservation metrics")
    parser.add_argument('--msa-dir', required=True, help='MSA directory to analyze')
    parser.add_argument('--output-file', required=True, help='Output CSV file')
    parser.add_argument('--strategy-name', default='Unknown', help='Strategy name for logging')
    
    args = parser.parse_args()
    
    # Analyze MSAs
    results_df = analyze_strategy_msas(args.msa_dir, args.strategy_name)
    
    if results_df.empty:
        print(f"No results generated for {args.strategy_name}")
        return
    
    # Save results
    os.makedirs(os.path.dirname(args.output_file), exist_ok=True)
    results_df.to_csv(args.output_file, index=False)
    
    print(f"\n✅ Enhanced summary saved to: {args.output_file}")
    print(f"📊 Analyzed {len(results_df)} genes")
    print(f"📈 Mean conservation: {results_df['mean_conservation'].mean():.3f}")
    print(f"📈 Mean identity: {results_df['mean_pairwise_identity'].mean():.3f}")

if __name__ == "__main__":
    main()
