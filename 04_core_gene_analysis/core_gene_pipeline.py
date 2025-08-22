#!/usr/bin/env python3
"""
Consolidated Core Gene Analysis Pipeline for E. faecalis genomes.
Combines all steps: identification, extraction, alignment, and diversity analysis.
"""

import os
import sys
import glob
import time
import math
import argparse
import subprocess
import pandas as pd
import numpy as np
from urllib.parse import unquote
from collections import Counter, defaultdict
from multiprocessing import Pool, cpu_count
from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import AlignInfo
import matplotlib.pyplot as plt
import re

# ============================================================================
# STEP 1: IDENTIFY CORE GENES
# ============================================================================

def parse_gff_attributes(attrs_string):
    """Parse GFF attributes with proper URL decoding."""
    attrs = {}
    for kv in attrs_string.split(';'):
        if '=' in kv:
            key, value = kv.split('=', 1)
            attrs[key] = unquote(value.strip())
    return attrs

def process_genome(prokka_dir):
    """Extract gene names and products from a single Prokka output."""
    genome_id = os.path.basename(prokka_dir)
    gff_file = os.path.join(prokka_dir, f"{genome_id}.gff")
    
    genes = set()
    gene_products = {}
    
    if os.path.exists(gff_file):
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 9:
                    attrs = parse_gff_attributes(parts[8])
                    
                    if parts[2] == 'CDS':
                        # Try gene first, then Name as fallback
                        gene_name = attrs.get('gene') or attrs.get('Name')
                        product = attrs.get('product', 'hypothetical protein')
                        
                        if gene_name:
                            gene_products[gene_name] = product
                            genes.add(gene_name)
                    
                    elif parts[2] == 'gene':
                        gene_name = attrs.get('gene') or attrs.get('Name')
                        if gene_name:
                            genes.add(gene_name)
    
    return genome_id, genes, gene_products

def identify_core_genes(prokka_dir, output_dir, threshold=0.95, threads=None):
    """Identify core genes present in ≥threshold of genomes."""
    
    print(f"\n{'='*60}")
    print("STEP 1: Identifying Core Genes")
    print(f"{'='*60}")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Find all Prokka directories
    prokka_dirs = glob.glob(os.path.join(prokka_dir, "ENT_*"))
    total_genomes = len(prokka_dirs)
    
    print(f"Found {total_genomes} Prokka output directories")
    print(f"Using {threads or cpu_count()} threads")
    
    # Process all genomes in parallel
    with Pool(threads or cpu_count()) as pool:
        results = pool.map(process_genome, prokka_dirs)
    
    # Aggregate results
    gene_counter = Counter()
    all_gene_products = {}
    genomes_processed = []
    
    for genome_id, genes, gene_products in results:
        if genes:
            genomes_processed.append(genome_id)
            for gene in genes:
                gene_counter[gene] += 1
                if gene in gene_products and gene not in all_gene_products:
                    all_gene_products[gene] = gene_products[gene]
    
    n_genomes = len(genomes_processed)
    print(f"Successfully processed {n_genomes} genomes")
    
    # Calculate prevalence and identify core genes
    prevalence_data = []
    core_genes = []
    min_count = math.ceil(n_genomes * threshold)  # Use ceiling to ensure proper threshold
    
    for gene, count in gene_counter.items():
        prevalence = count / n_genomes
        prevalence_data.append({
            'gene': gene,
            'count': count,
            'prevalence': prevalence,
            'product': all_gene_products.get(gene, 'Unknown')
        })
        
        if count >= min_count:
            core_genes.append(gene)
    
    # Save results
    df = pd.DataFrame(prevalence_data)
    df = df.sort_values('prevalence', ascending=False)
    df.to_csv(os.path.join(output_dir, 'gene_prevalence_stats.csv'), index=False)
    
    # Save core genes list
    core_genes_file = os.path.join(output_dir, f'core_genes_{int(threshold*100)}pct.txt')
    with open(core_genes_file, 'w') as f:
        for gene in sorted(core_genes):
            f.write(f"{gene}\n")
    
    print(f"\nResults:")
    print(f"  Total genomes found: {total_genomes}")
    print(f"  Genomes successfully processed: {n_genomes}")
    if total_genomes != n_genomes:
        print(f"  ⚠️  Warning: {total_genomes - n_genomes} genomes failed to process")
    print(f"  Total unique genes: {len(gene_counter)}")
    print(f"  Core genes (≥{threshold*100}% of {n_genomes} genomes): {len(core_genes)}")
    print(f"  Minimum gene count for core: {min_count}")
    print(f"  Saved to: {core_genes_file}")
    
    return core_genes, n_genomes

# ============================================================================
# STEP 2: EXTRACT SEQUENCES
# ============================================================================

def extract_gene_from_genome(args):
    """Extract a specific gene from a genome's Prokka output."""
    gene_name, genome_dir = args
    
    genome_id = os.path.basename(genome_dir.rstrip('/'))
    
    gff_file = os.path.join(genome_dir, f"{genome_id}.gff")
    fna_file = os.path.join(genome_dir, f"{genome_id}.fna")
    
    if not os.path.exists(gff_file) or not os.path.exists(fna_file):
        return None
    
    # Parse genome sequence
    genome_seqs = {}
    try:
        for record in SeqIO.parse(fna_file, "fasta"):
            genome_seqs[record.id] = record
    except:
        return None
    
    # Find gene coordinates in GFF
    gene_coords = None
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 9 and parts[2] == 'gene':
                attrs = parse_gff_attributes(parts[8])
                # Try gene first, then Name as fallback
                found_gene = attrs.get('gene') or attrs.get('Name')
                if found_gene == gene_name:
                    gene_coords = {
                        'contig': parts[0],
                        'start': int(parts[3]),
                        'end': int(parts[4]),
                        'strand': parts[6]
                    }
                    break
    
    if not gene_coords or gene_coords['contig'] not in genome_seqs:
        return None
    
    # Extract sequence
    contig_seq = genome_seqs[gene_coords['contig']].seq
    gene_seq = contig_seq[gene_coords['start']-1:gene_coords['end']]
    
    if gene_coords['strand'] == '-':
        gene_seq = gene_seq.reverse_complement()
    
    # Create SeqRecord
    record = SeqRecord(
        gene_seq,
        id=genome_id,
        description=f"{gene_name} from {genome_id}"
    )
    
    return record

def extract_sequences(core_genes, prokka_dir, output_dir, threads=None):
    """Extract sequences for all core genes from all genomes."""
    
    print(f"\n{'='*60}")
    print("STEP 2: Extracting Core Gene Sequences")
    print(f"{'='*60}")
    
    seq_dir = os.path.join(output_dir, 'core_gene_sequences')
    os.makedirs(seq_dir, exist_ok=True)
    
    prokka_dirs = glob.glob(os.path.join(prokka_dir, "ENT_*"))
    print(f"Processing {len(core_genes)} genes from {len(prokka_dirs)} genomes")
    
    success_count = 0
    
    for i, gene in enumerate(sorted(core_genes), 1):
        if i % 100 == 0:
            print(f"  Processing gene {i}/{len(core_genes)}...")
        
        # Create arguments for parallel processing
        args_list = [(gene, genome_dir) for genome_dir in prokka_dirs]
        
        # Extract sequences in parallel
        with Pool(threads or cpu_count()) as pool:
            sequences = pool.map(extract_gene_from_genome, args_list)
        
        # Filter out None results
        sequences = [seq for seq in sequences if seq is not None]
        
        if sequences:
            output_file = os.path.join(seq_dir, f"{gene}.fasta")
            SeqIO.write(sequences, output_file, "fasta")
            success_count += 1
    
    print(f"\nExtracted sequences for {success_count}/{len(core_genes)} genes")
    
    # Create extraction summary file
    summary_file = os.path.join(output_dir, 'extraction_summary.txt')
    with open(summary_file, 'w') as f:
        f.write("Core Gene Sequence Extraction Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total core genes: {len(core_genes)}\n")
        f.write(f"Successfully extracted: {success_count}\n")
        f.write(f"Failed to extract: {len(core_genes) - success_count}\n")
        f.write(f"Genomes processed: {len(prokka_dirs)}\n")
        f.write(f"\nOutput directory: {seq_dir}\n")
        
        if success_count < len(core_genes):
            f.write(f"\nGenes that failed extraction:\n")
            extracted_genes = set(f.replace('.fasta', '') for f in os.listdir(seq_dir) if f.endswith('.fasta'))
            failed_genes = set(core_genes) - extracted_genes
            for gene in sorted(failed_genes):
                f.write(f"  - {gene}\n")
    
    print(f"  Saved extraction summary to: {summary_file}")
    
    return seq_dir

# ============================================================================
# STEP 3: CREATE MULTIPLE SEQUENCE ALIGNMENTS
# ============================================================================

def create_msa_for_gene(args):
    """Create MSA for a single gene using MAFFT."""
    gene_file, output_dir, timeout = args
    
    gene_name = os.path.basename(gene_file).replace('.fasta', '')
    output_file = os.path.join(output_dir, f"{gene_name}_aligned.fasta")
    
    # Check if input file has sequences
    try:
        sequences = list(SeqIO.parse(gene_file, "fasta"))
        if len(sequences) < 2:
            return gene_name, f"Skipped (only {len(sequences)} sequences)"
    except:
        return gene_name, "Error reading sequences"
    
    # Run MAFFT
    try:
        cmd = ["mafft", "--auto", "--quiet", gene_file]
        actual_timeout = None if timeout == 0 else timeout
        with open(output_file, 'w') as outf:
            result = subprocess.run(cmd, stdout=outf, stderr=subprocess.PIPE, 
                                    text=True, timeout=actual_timeout)
        
        if result.returncode == 0:
            try:
                alignment = AlignIO.read(output_file, "fasta")
                return gene_name, f"Success ({len(alignment)} sequences, {alignment.get_alignment_length()} bp)"
            except:
                return gene_name, "Error: Invalid alignment produced"
        else:
            return gene_name, f"MAFFT failed: {result.stderr}"
            
    except subprocess.TimeoutExpired:
        return gene_name, f"Timeout (>{timeout} seconds)"
    except FileNotFoundError:
        return gene_name, "Error: MAFFT not found"
    except Exception as e:
        return gene_name, f"Error: {str(e)}"

def create_alignments(seq_dir, output_dir, threads=None, timeout=0):
    """Create MSAs for all core genes."""
    
    print(f"\n{'='*60}")
    print("STEP 3: Creating Multiple Sequence Alignments")
    print(f"{'='*60}")
    
    align_dir = os.path.join(output_dir, 'core_gene_alignments')
    os.makedirs(align_dir, exist_ok=True)
    
    # Find all sequence files
    fasta_files = glob.glob(os.path.join(seq_dir, "*.fasta"))
    print(f"Creating alignments for {len(fasta_files)} genes")
    print(f"Using {threads or cpu_count()} threads")
    print(f"Timeout per gene: {timeout} seconds")
    
    # Create arguments for parallel processing
    args_list = [(f, align_dir, timeout) for f in fasta_files]
    
    # Process in parallel
    success_count = 0
    failed_genes = []
    
    with Pool(threads or cpu_count()) as pool:
        results = pool.map(create_msa_for_gene, args_list)
    
    # Process results
    for gene_name, status in results:
        if "Success" in status:
            success_count += 1
        else:
            failed_genes.append((gene_name, status))
    
    print(f"\nAlignment Results:")
    print(f"  Successful: {success_count}/{len(fasta_files)}")
    
    if failed_genes:
        print(f"  Failed: {len(failed_genes)}")
        for gene, reason in failed_genes[:5]:
            print(f"    - {gene}: {reason}")
        if len(failed_genes) > 5:
            print(f"    ... and {len(failed_genes)-5} more")
    
    # Save summary
    summary_file = os.path.join(align_dir, 'msa_summary.txt')
    with open(summary_file, 'w') as f:
        f.write(f"MSA Creation Summary\n")
        f.write(f"====================\n")
        f.write(f"Total genes: {len(fasta_files)}\n")
        f.write(f"Successful: {success_count}\n")
        f.write(f"Failed: {len(failed_genes)}\n\n")
        
        if failed_genes:
            f.write("Failed alignments:\n")
            for gene, reason in failed_genes:
                f.write(f"  {gene}: {reason}\n")
    
    return align_dir

# ============================================================================
# STEP 4: CALCULATE DIVERSITY METRICS
# ============================================================================

def calculate_conservation_metrics(alignment_file):
    """Calculate conservation metrics for a single gene alignment."""
    
    gene_name = os.path.basename(alignment_file).replace('_aligned.fasta', '')
    
    try:
        alignment = AlignIO.read(alignment_file, "fasta")
        n_sequences = len(alignment)
        alignment_length = alignment.get_alignment_length()
        
        if n_sequences < 2:
            return {
                'gene': gene_name,
                'error': 'Too few sequences',
                'n_sequences': n_sequences
            }
        
        # Calculate position-wise conservation scores
        position_scores = []
        gap_counts = []
        
        for pos in range(alignment_length):
            column = alignment[:, pos]
            
            # Count characters (excluding gaps)
            char_counts = Counter([c for c in column if c != '-'])
            gap_count = column.count('-')
            gap_counts.append(gap_count)
            
            if len(char_counts) == 0:
                conservation_score = 0
            else:
                # Calculate Shannon entropy
                total = sum(char_counts.values())
                entropy = -sum((count/total) * np.log2(count/total) 
                             for count in char_counts.values() if count > 0)
                
                # Convert entropy to conservation (max entropy for DNA = 2)
                # Clamp to [0, 1] to handle edge cases with ambiguous bases
                max_entropy = 2.0
                conservation_score = 1 - (entropy / max_entropy)
                conservation_score = max(0.0, min(1.0, conservation_score))
            
            position_scores.append(conservation_score)
        
        # Calculate overall metrics
        mean_conservation = np.mean(position_scores)
        median_conservation = np.median(position_scores)
        std_conservation = np.std(position_scores)
        
        # Calculate percentage of highly conserved positions
        highly_conserved = sum(1 for s in position_scores if s > 0.9) / len(position_scores)
        moderately_conserved = sum(1 for s in position_scores if 0.5 <= s <= 0.9) / len(position_scores)
        variable = sum(1 for s in position_scores if s < 0.5) / len(position_scores)
        
        # Gap statistics
        mean_gaps = np.mean(gap_counts)
        max_gaps = max(gap_counts)
        positions_with_gaps = sum(1 for g in gap_counts if g > 0) / len(gap_counts)
        
        return {
            'gene': gene_name,
            'n_sequences': n_sequences,
            'alignment_length': alignment_length,
            'mean_conservation': mean_conservation,
            'median_conservation': median_conservation,
            'std_conservation': std_conservation,
            'highly_conserved_pct': highly_conserved * 100,
            'moderately_conserved_pct': moderately_conserved * 100,
            'variable_pct': variable * 100,
            'mean_gaps_per_position': mean_gaps,
            'max_gaps_per_position': max_gaps,
            'positions_with_gaps_pct': positions_with_gaps * 100
        }
        
    except Exception as e:
        return {
            'gene': gene_name,
            'error': str(e)
        }

def calculate_diversity(align_dir, output_dir, threads=None):
    """Calculate diversity metrics for all aligned genes."""
    
    print(f"\n{'='*60}")
    print("STEP 4: Calculating Diversity Metrics")
    print(f"{'='*60}")
    
    # Find all alignment files
    alignment_files = glob.glob(os.path.join(align_dir, "*_aligned.fasta"))
    print(f"Processing {len(alignment_files)} alignments")
    
    # Process in parallel
    with Pool(threads or cpu_count()) as pool:
        results = pool.map(calculate_conservation_metrics, alignment_files)
    
    # Convert to DataFrame
    df = pd.DataFrame(results)
    
    # Separate successful and failed
    # Check if 'error' column exists
    if 'error' in df.columns:
        successful = df[~df['error'].notna()]
        failed = df[df['error'].notna()]
    else:
        # All successful - no error column exists
        successful = df
        failed = pd.DataFrame()
    
    print(f"\nProcessed {len(successful)} alignments successfully")
    if len(failed) > 0:
        print(f"Failed to process {len(failed)} alignments")
    
    # Save results
    output_file = os.path.join(output_dir, 'core_gene_conservation_metrics.csv')
    if len(successful) > 0:
        successful.to_csv(output_file, index=False)
    else:
        print("Warning: No successful alignments to save")
        return pd.DataFrame()
    
    # Create summary statistics
    summary_file = os.path.join(output_dir, 'core_gene_conservation_summary.txt')
    with open(summary_file, 'w') as f:
        f.write("Core Gene Conservation Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total genes analyzed: {len(successful)}\n")
        f.write(f"Average sequences per gene: {successful['n_sequences'].mean():.1f}\n")
        f.write(f"Average alignment length: {successful['alignment_length'].mean():.1f} bp\n\n")
        
        f.write("Conservation Statistics:\n")
        f.write(f"  Mean conservation: {successful['mean_conservation'].mean():.3f} ± {successful['mean_conservation'].std():.3f}\n")
        f.write(f"  Highly conserved positions: {successful['highly_conserved_pct'].mean():.1f}%\n")
        f.write(f"  Moderately conserved: {successful['moderately_conserved_pct'].mean():.1f}%\n")
        f.write(f"  Variable positions: {successful['variable_pct'].mean():.1f}%\n\n")
        
        f.write("Gap Statistics:\n")
        f.write(f"  Average gaps per position: {successful['mean_gaps_per_position'].mean():.2f}\n")
        f.write(f"  Positions with gaps: {successful['positions_with_gaps_pct'].mean():.1f}%\n")
    
    print(f"\nResults saved to:")
    print(f"  - {output_file}")
    print(f"  - {summary_file}")
    
    return successful

# ============================================================================
# STEP 5: THRESHOLD ANALYSIS
# ============================================================================

def analyze_thresholds(prokka_dir, output_dir):
    """Analyze how gene count changes with different prevalence thresholds."""
    
    print(f"\n{'='*60}")
    print("STEP 5: Threshold Analysis")
    print(f"{'='*60}")
    
    # Load prevalence data
    prevalence_file = os.path.join(output_dir, 'gene_prevalence_stats.csv')
    if not os.path.exists(prevalence_file):
        print("Error: gene_prevalence_stats.csv not found. Run step 1 first.")
        return
    
    df = pd.read_csv(prevalence_file)
    
    # Calculate threshold curve
    thresholds = np.arange(0.0, 1.01, 0.01)
    gene_counts = []
    
    for threshold in thresholds:
        count = (df['prevalence'] >= threshold).sum()
        gene_counts.append(count)
    
    # Create plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Full range plot
    ax1.plot(thresholds * 100, gene_counts, 'b-', linewidth=2)
    ax1.axvline(x=95, color='r', linestyle='--', label='95% threshold')
    ax1.axhline(y=gene_counts[95], color='r', linestyle='--', alpha=0.5)
    ax1.set_xlabel('Prevalence Threshold (%)', fontsize=12)
    ax1.set_ylabel('Number of Core Genes', fontsize=12)
    ax1.set_title('Core Gene Count vs Prevalence Threshold', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Zoomed plot (90-100%)
    high_thresholds = thresholds[90:]
    high_counts = gene_counts[90:]
    
    ax2.plot(high_thresholds * 100, high_counts, 'b-', linewidth=2)
    ax2.axvline(x=95, color='r', linestyle='--', label='95% threshold')
    ax2.axhline(y=gene_counts[95], color='r', linestyle='--', alpha=0.5)
    ax2.set_xlabel('Prevalence Threshold (%)', fontsize=12)
    ax2.set_ylabel('Number of Core Genes', fontsize=12)
    ax2.set_title('Core Gene Count (90-100% range)', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    # Add annotations
    ax2.annotate(f'{gene_counts[95]} genes at 95%', 
                xy=(95, gene_counts[95]), 
                xytext=(96, gene_counts[95] + 50),
                arrowprops=dict(arrowstyle='->', color='red'),
                fontsize=10)
    
    plt.tight_layout()
    
    # Save plots
    plot_file = os.path.join(output_dir, 'core_gene_threshold_curve.png')
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"Saved threshold analysis plot to: {plot_file}")
    
    # Save threshold data
    threshold_df = pd.DataFrame({
        'threshold_pct': thresholds * 100,
        'gene_count': gene_counts
    })
    threshold_file = os.path.join(output_dir, 'core_gene_threshold_summary.csv')
    threshold_df.to_csv(threshold_file, index=False)
    
    # Print key statistics
    print(f"\nThreshold Analysis Results:")
    print(f"  100% prevalence: {gene_counts[100]} genes")
    print(f"  99% prevalence: {gene_counts[99]} genes")
    print(f"  95% prevalence: {gene_counts[95]} genes")
    print(f"  90% prevalence: {gene_counts[90]} genes")
    print(f"  50% prevalence: {gene_counts[50]} genes")

# ============================================================================
# MAIN PIPELINE
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Core Gene Analysis Pipeline for E. faecalis genomes',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Pipeline Steps:
  1. Identify core genes (≥95% prevalence)
  2. Extract sequences from Prokka output
  3. Create MSAs using MAFFT
  4. Calculate diversity metrics
  5. Threshold analysis and plots

Examples:
  %(prog)s --prokka-dir ../01_prokka_annotation/prokka_output
  %(prog)s --start-step 3 --threads 40
  %(prog)s --start-step 5  # Only run threshold analysis
        """
    )
    
    parser.add_argument('--prokka-dir', 
                       default='../01_prokka_annotation/output/prokka_results',
                       help='Directory containing Prokka output folders')
    parser.add_argument('--output-dir', 
                       default='output',
                       help='Output directory (default: output)')
    parser.add_argument('--threshold', 
                       type=float, 
                       default=0.95,
                       help='Prevalence threshold for core genes (default: 0.95)')
    parser.add_argument('--threads', 
                       type=int,
                       help='Number of threads (default: all available)')
    parser.add_argument('--start-step', 
                       type=int, 
                       default=1,
                       choices=[1, 2, 3, 4, 5],
                       help='Start from specific step (1-5)')
    parser.add_argument('--mafft-timeout', 
                       type=int, 
                       default=0,
                       help='Timeout for MAFFT per gene in seconds (default: 0=no timeout)')
    
    args = parser.parse_args()
    
    # Set environment for MAFFT
    if 'MAFFT_TMPDIR' not in os.environ:
        os.environ['MAFFT_TMPDIR'] = '/vol/tmp'
    
    print(f"\n{'='*60}")
    print("CORE GENE ANALYSIS PIPELINE")
    print(f"{'='*60}")
    print(f"Prokka directory: {args.prokka_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Threshold: {args.threshold*100}%")
    print(f"Threads: {args.threads or cpu_count()}")
    print(f"Starting from step: {args.start_step}")
    
    # Check prokka directory exists
    if not os.path.exists(args.prokka_dir):
        print(f"\nError: Prokka directory not found: {args.prokka_dir}")
        sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Run pipeline steps
    core_genes = None
    seq_dir = None
    align_dir = None
    
    # Step 1: Identify core genes
    if args.start_step <= 1:
        core_genes, n_genomes = identify_core_genes(
            args.prokka_dir, 
            args.output_dir, 
            args.threshold,
            args.threads
        )
    else:
        # Load existing core genes
        core_genes_file = os.path.join(args.output_dir, f'core_genes_{int(args.threshold*100)}pct.txt')
        if os.path.exists(core_genes_file):
            with open(core_genes_file, 'r') as f:
                core_genes = [line.strip() for line in f if line.strip()]
            print(f"\nLoaded {len(core_genes)} existing core genes")
        else:
            print(f"\nError: Core genes file not found. Run from step 1.")
            sys.exit(1)
    
    # Step 2: Extract sequences
    if args.start_step <= 2:
        if core_genes:
            seq_dir = extract_sequences(
                core_genes,
                args.prokka_dir,
                args.output_dir,
                args.threads
            )
    else:
        seq_dir = os.path.join(args.output_dir, 'core_gene_sequences')
        if os.path.exists(seq_dir):
            print(f"\nUsing existing sequences in: {seq_dir}")
        else:
            print(f"\nError: Sequence directory not found. Run from step 2.")
            sys.exit(1)
    
    # Step 3: Create MSAs
    if args.start_step <= 3:
        if seq_dir:
            align_dir = create_alignments(
                seq_dir,
                args.output_dir,
                args.threads,
                args.mafft_timeout
            )
    else:
        align_dir = os.path.join(args.output_dir, 'core_gene_alignments')
        if os.path.exists(align_dir):
            print(f"\nUsing existing alignments in: {align_dir}")
        else:
            print(f"\nError: Alignment directory not found. Run from step 3.")
            sys.exit(1)
    
    # Step 4: Calculate diversity
    if args.start_step <= 4:
        if align_dir:
            diversity_df = calculate_diversity(
                align_dir,
                args.output_dir,
                args.threads
            )
    
    # Step 5: Threshold analysis
    if args.start_step <= 5:
        analyze_thresholds(args.prokka_dir, args.output_dir)
    
    print(f"\n{'='*60}")
    print("PIPELINE COMPLETE!")
    print(f"{'='*60}")
    print(f"Results saved to: {args.output_dir}/")
    print("\nKey output files:")
    print(f"  - gene_prevalence_stats.csv")
    print(f"  - core_genes_{int(args.threshold*100)}pct.txt")
    print(f"  - extraction_summary.txt")
    print(f"  - core_gene_sequences/")
    print(f"  - core_gene_alignments/")
    print(f"  - core_gene_conservation_metrics.csv")
    print(f"  - core_gene_threshold_curve.png")

if __name__ == "__main__":
    main()