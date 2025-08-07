#!/usr/bin/env python3
"""
Calculate conservation scores and diversity metrics for core genes.
Creates comprehensive diversity analysis from MSA files.
"""

import os
import glob
import pandas as pd
import numpy as np
from Bio import AlignIO
from Bio.Align import AlignInfo
from collections import defaultdict, Counter
import argparse
import time
from multiprocessing import Pool, cpu_count

def calculate_conservation_metrics(alignment_file):
    """Calculate various conservation metrics for a single gene alignment."""
    
    gene_name = os.path.basename(alignment_file).replace('_aligned.fasta', '')
    
    try:
        # Read alignment
        alignment = AlignIO.read(alignment_file, "fasta")
        n_sequences = len(alignment)
        alignment_length = alignment.get_alignment_length()
        
        if n_sequences < 2:
            return {
                'gene': gene_name,
                'error': 'Too few sequences',
                'n_sequences': n_sequences
            }
        
        # Calculate position-wise conservation scores (Shannon entropy)
        position_scores = []
        gap_counts = []
        
        for pos in range(alignment_length):
            column = alignment[:, pos]
            
            # Count characters (excluding gaps for conservation)
            char_counts = Counter([c for c in column if c != '-'])
            gap_count = column.count('-')
            gap_counts.append(gap_count)
            
            if len(char_counts) == 0:
                # All gaps
                conservation_score = 0
            else:
                # Shannon entropy based conservation
                total_chars = sum(char_counts.values())
                entropy = 0
                for count in char_counts.values():
                    if count > 0:
                        p = count / total_chars
                        entropy -= p * np.log2(p)
                
                # Convert to conservation score (0-1, higher = more conserved)
                max_entropy = np.log2(4)  # For DNA (A,T,G,C)
                conservation_score = 1 - (entropy / max_entropy) if max_entropy > 0 else 0
            
            position_scores.append(conservation_score)
        
        # Overall metrics
        mean_conservation = np.mean(position_scores)
        std_conservation = np.std(position_scores)
        min_conservation = np.min(position_scores)
        max_conservation = np.max(position_scores)
        
        # Gap statistics
        mean_gaps = np.mean(gap_counts)
        total_gaps = sum(gap_counts)
        gap_percentage = (total_gaps / (n_sequences * alignment_length)) * 100
        
        # Sequence identity statistics
        identities = []
        for i in range(n_sequences):
            for j in range(i + 1, n_sequences):
                seq1 = str(alignment[i].seq).replace('-', '')
                seq2 = str(alignment[j].seq).replace('-', '')
                
                if len(seq1) > 0 and len(seq2) > 0:
                    # Simple identity calculation
                    min_len = min(len(seq1), len(seq2))
                    matches = sum(1 for k in range(min_len) if seq1[k] == seq2[k])
                    identity = matches / min_len if min_len > 0 else 0
                    identities.append(identity)
        
        mean_identity = np.mean(identities) if identities else 0
        std_identity = np.std(identities) if identities else 0
        min_identity = np.min(identities) if identities else 0
        max_identity = np.max(identities) if identities else 0
        
        # Most conserved regions (windows of 10 positions)
        window_size = 10
        window_scores = []
        for i in range(0, len(position_scores) - window_size + 1, window_size):
            window_score = np.mean(position_scores[i:i + window_size])
            window_scores.append(window_score)
        
        most_conserved_window = np.max(window_scores) if window_scores else 0
        least_conserved_window = np.min(window_scores) if window_scores else 0
        
        return {
            'gene': gene_name,
            'n_sequences': n_sequences,
            'alignment_length': alignment_length,
            'mean_conservation': mean_conservation,
            'std_conservation': std_conservation,
            'min_conservation': min_conservation,
            'max_conservation': max_conservation,
            'mean_pairwise_identity': mean_identity,
            'std_pairwise_identity': std_identity,
            'min_pairwise_identity': min_identity,
            'max_pairwise_identity': max_identity,
            'gap_percentage': gap_percentage,
            'mean_gaps_per_position': mean_gaps,
            'most_conserved_window': most_conserved_window,
            'least_conserved_window': least_conserved_window,
            'conservation_category': categorize_conservation(mean_conservation),
            'error': None
        }
        
    except Exception as e:
        return {
            'gene': gene_name,
            'error': str(e),
            'n_sequences': 0
        }

def categorize_conservation(score):
    """Categorize conservation score into qualitative bins."""
    if score >= 0.9:
        return 'Highly conserved'
    elif score >= 0.7:
        return 'Moderately conserved'
    elif score >= 0.5:
        return 'Lowly conserved'
    else:
        return 'Highly variable'

def process_gene_alignment(alignment_file):
    """Wrapper for multiprocessing."""
    return calculate_conservation_metrics(alignment_file)

def main():
    parser = argparse.ArgumentParser(description="Calculate conservation metrics for core genes")
    parser.add_argument("--alignments-dir", default="output/core_gene_alignments",
                        help="Directory containing gene alignments")
    parser.add_argument("--output-dir", default="output",
                        help="Output directory")
    parser.add_argument("--threads", type=int, default=None,
                        help="Number of threads to use (default: use SLURM_CPUS_PER_TASK if set, else all CPUs)")
    parser.add_argument("--progress-every", type=int, default=50,
                        help="Print progress every N genes (default: 50)")
    parser.add_argument("--checkpoint-every", type=int, default=0,
                        help="Write partial results every N genes (0 = disabled)")
    
    args = parser.parse_args()
    
    # Check input directory
    if not os.path.exists(args.alignments_dir):
        print(f"Error: Alignments directory not found: {args.alignments_dir}")
        print("Run create_core_gene_msa.py first!")
        return 1
    
    # Find alignment files
    alignment_files = glob.glob(os.path.join(args.alignments_dir, "*_aligned.fasta"))
    
    if not alignment_files:
        print(f"Error: No alignment files found in {args.alignments_dir}")
        return 1
    
    print(f"{'='*60}")
    print(f"Core Gene Conservation Analysis")
    print(f"{'='*60}")
    total_genes = len(alignment_files)
    print(f"Alignments directory: {args.alignments_dir}")
    print(f"Number of genes: {total_genes}")
    # Resolve threads: prefer CLI, then SLURM_CPUS_PER_TASK, else all CPUs
    if args.threads is None:
        try:
            from os import environ
            slurm_threads = int(environ.get("SLURM_CPUS_PER_TASK", "0"))
        except ValueError:
            slurm_threads = 0
        args.threads = slurm_threads if slurm_threads > 0 else cpu_count()

    print(f"Threads: {args.threads}")
    print(f"{'='*60}\n")
    
    # Process alignments in parallel with progress logging
    print("Calculating conservation metrics...", flush=True)
    start_time = time.time()
    results = []
    failed_results = []
    successful_results = []
    processed = 0
    next_report = args.progress_every
    partial_file = None
    if args.checkpoint_every and args.checkpoint_every > 0:
        os.makedirs(args.output_dir, exist_ok=True)
        partial_file = os.path.join(args.output_dir, 'core_gene_conservation_metrics.partial.csv')

    # Use imap_unordered for smoother progress
    with Pool(args.threads) as pool:
        for res in pool.imap_unordered(process_gene_alignment, alignment_files, chunksize=5):
            results.append(res)
            processed += 1
            if res.get('error') is None:
                successful_results.append(res)
            else:
                failed_results.append(res)

            # Periodic progress report
            if processed >= next_report or processed == total_genes:
                elapsed = time.time() - start_time
                rate = processed / elapsed if elapsed > 0 else 0
                remaining = total_genes - processed
                eta_sec = remaining / rate if rate > 0 else 0
                eta_min = int(eta_sec // 60)
                eta_sec_rem = int(eta_sec % 60)
                print(f"Progress: {processed}/{total_genes} (" \
                      f"{processed/total_genes*100:.1f}%) | " \
                      f"Elapsed: {int(elapsed//60)}m{int(elapsed%60)}s | " \
                      f"ETA: {eta_min}m{eta_sec_rem}s | " \
                      f"OK: {len(successful_results)} | Fail: {len(failed_results)}",
                      flush=True)
                next_report += args.progress_every

                # Optional checkpoint write
                if partial_file and (processed % args.checkpoint_every == 0 or processed == total_genes):
                    try:
                        pd.DataFrame([r for r in successful_results]).to_csv(partial_file, index=False)
                        print(f"Checkpoint written: {partial_file} (rows: {len(successful_results)})", flush=True)
                    except Exception as e:
                        print(f"Warning: failed to write checkpoint: {e}", flush=True)
    
    # Results already split during progress loop

    if failed_results:
        print(f"\n⚠️  Warning: {len(failed_results)} genes failed analysis:")
        for result in failed_results[:5]:  # Show first 5 errors
            print(f"  {result['gene']}: {result['error']}")
        if len(failed_results) > 5:
            print(f"  ... and {len(failed_results) - 5} more")
    
    if not successful_results:
        print("❌ No genes successfully analyzed")
        return 1
    
    # Create DataFrame
    df = pd.DataFrame(successful_results)
    
    # Sort by conservation score
    df = df.sort_values('mean_conservation', ascending=False)
    
    # Save results
    os.makedirs(args.output_dir, exist_ok=True)
    output_file = os.path.join(args.output_dir, 'core_gene_conservation_metrics.csv')
    df.to_csv(output_file, index=False, float_format='%.6f')
    
    # Create summary statistics
    summary_stats = {
        'total_genes_analyzed': len(successful_results),
        'mean_conservation_score': df['mean_conservation'].mean(),
        'std_conservation_score': df['mean_conservation'].std(),
        'highly_conserved_genes': (df['conservation_category'] == 'Highly conserved').sum(),
        'moderately_conserved_genes': (df['conservation_category'] == 'Moderately conserved').sum(),
        'lowly_conserved_genes': (df['conservation_category'] == 'Lowly conserved').sum(),
        'highly_variable_genes': (df['conservation_category'] == 'Highly variable').sum(),
        'mean_pairwise_identity': df['mean_pairwise_identity'].mean(),
        'mean_gap_percentage': df['gap_percentage'].mean(),
    }
    
    # Save summary
    summary_file = os.path.join(args.output_dir, 'core_gene_conservation_summary.txt')
    with open(summary_file, 'w') as f:
        f.write("Core Gene Conservation Analysis Summary\n")
        f.write("======================================\n\n")
        for key, value in summary_stats.items():
            if isinstance(value, float):
                f.write(f"{key}: {value:.4f}\n")
            else:
                f.write(f"{key}: {value}\n")
        
        f.write(f"\nTop 10 most conserved genes:\n")
        for _, row in df.head(10).iterrows():
            f.write(f"  {row['gene']}: {row['mean_conservation']:.4f} ({row['conservation_category']})\n")
        
        f.write(f"\nTop 10 most variable genes:\n")
        for _, row in df.tail(10).iterrows():
            f.write(f"  {row['gene']}: {row['mean_conservation']:.4f} ({row['conservation_category']})\n")
    
    # Display results
    print(f"\n{'='*60}")
    print("Conservation Analysis Results:")
    print(f"{'='*60}")
    print(f"Successfully analyzed: {len(successful_results)} genes")
    print(f"Failed: {len(failed_results)} genes")
    print(f"\nConservation distribution:")
    print(f"  Highly conserved (≥0.9):   {summary_stats['highly_conserved_genes']} genes")
    print(f"  Moderately conserved (≥0.7): {summary_stats['moderately_conserved_genes']} genes")
    print(f"  Lowly conserved (≥0.5):     {summary_stats['lowly_conserved_genes']} genes")
    print(f"  Highly variable (<0.5):     {summary_stats['highly_variable_genes']} genes")
    print(f"\nOverall statistics:")
    print(f"  Mean conservation score: {summary_stats['mean_conservation_score']:.4f}")
    print(f"  Mean pairwise identity:  {summary_stats['mean_pairwise_identity']:.4f}")
    print(f"  Mean gap percentage:     {summary_stats['mean_gap_percentage']:.2f}%")
    
    print(f"\n📁 Results saved to:")
    print(f"  {output_file}")
    print(f"  {summary_file}")
    
    print(f"\n✅ Conservation analysis complete!")
    return 0

if __name__ == "__main__":
    exit(main())