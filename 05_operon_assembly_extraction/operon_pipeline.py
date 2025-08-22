#!/usr/bin/env python3
"""
Consolidated operon extraction and analysis pipeline.
Combines all functionality from separate scripts into a single, maintainable module.
"""

import os
import sys
import time
import gzip
import re
import argparse
import subprocess
import multiprocessing
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for headless plotting
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Optional imports
try:
    import logomaker
    LOGOMAKER_AVAILABLE = True
except ImportError:
    LOGOMAKER_AVAILABLE = False
    print("Warning: logomaker not available, sequence logos will be skipped")


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def _normalize_variants(s: str) -> set:
    """Generate multiple ID variants for robust contig matching."""
    raw = s.split()[0]
    no_prefix = raw.split('|', 2)[2] if raw.startswith('gnl|') and raw.count('|') >= 2 else raw
    no_ver_raw = re.sub(r'\.\d+$', '', raw)
    no_ver_no_prefix = re.sub(r'\.\d+$', '', no_prefix)
    return {raw, no_prefix, no_ver_raw, no_ver_no_prefix}


def _infer_genome_id_from_blast_filename(filename: str) -> str:
    """Infer genome ID from BLAST result filename."""
    base = os.path.basename(filename)
    
    # Remove BLAST suffixes
    for suf in ("_genes_blast.txt", "_noncoding_blast.txt", "_blast.txt", ".blast", ".txt"):
        if base.endswith(suf):
            base = base[:-len(suf)]
            break
    
    # Remove .fasta if present (for cases like ENT_CA1914AA_AS.genome.fasta)
    if base.endswith(".fasta"):
        base = base[:-6]  # Remove ".fasta"
    
    return base


def parse_blast_hit(blast_line: str) -> dict:
    """Parse a BLAST hit line."""
    fields = blast_line.strip().split('\t')
    return {
        'qseqid': fields[0],
        'sseqid': fields[1],
        'pident': float(fields[2]),
        'length': int(fields[3]),
        'sstart': int(fields[8]),
        'send': int(fields[9]),
        'qcovs': float(fields[12]) if len(fields) > 12 else 100.0
    }


# ============================================================================
# SEQUENCE EXTRACTION MODULE
# ============================================================================

def extract_operon_sequences(blast_dir: str, output_dir: str, 
                            prokka_dir: Optional[str] = None,
                            assemblies_dir: Optional[str] = None,
                            min_identity: float = 90.0,
                            min_coverage: float = 80.0,
                            source: str = 'assemblies',
                            require_complete: bool = False) -> None:
    """
    Extract operon gene sequences from genomes based on BLAST results.
    
    Args:
        blast_dir: Directory containing BLAST result files
        output_dir: Directory to write extracted sequences
        prokka_dir: Directory with Prokka annotations (for source='prokka')
        assemblies_dir: Directory with genome assemblies (for source='assemblies')
        min_identity: Minimum sequence identity threshold
        min_coverage: Minimum query coverage threshold
        source: 'prokka' or 'assemblies'
        require_complete: Only extract from genomes with all 7 operon genes
    """
    print(f"\n{'='*60}")
    print(f"Extracting operon sequences (source={source})")
    print(f"{'='*60}")
    
    # Operon genes
    operon_genes = ['frpC', 'glpC', 'ptsD', 'ptsC', 'ptsB', 'ptsA', 'fruR']
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Process BLAST results
    print(f"Reading BLAST results from {blast_dir}")
    # Look specifically for gene BLAST files
    blast_files = [f for f in os.listdir(blast_dir) if "_genes_blast.txt" in f]
    
    if not blast_files:
        print(f"Error: No BLAST files found in {blast_dir}")
        sys.exit(2)
    
    print(f"Found {len(blast_files)} BLAST files to process")
    
    genome_best_hits = defaultdict(dict)
    
    for blast_file in blast_files:
        genome_id = _infer_genome_id_from_blast_filename(blast_file)
        file_path = os.path.join(blast_dir, blast_file)
        
        with open(file_path, 'r') as f:
            for line in f:
                if not line.strip():
                    continue
                
                hit = parse_blast_hit(line)
                # Extract gene name from format like operon_CDS_1|frpC|product
                parts = hit['qseqid'].split('|')
                gene = parts[1] if len(parts) > 1 else parts[0]
                
                if gene not in operon_genes:
                    continue
                
                if hit['pident'] < min_identity or hit['qcovs'] < min_coverage:
                    continue
                
                # Keep best hit per gene
                if gene not in genome_best_hits[genome_id]:
                    genome_best_hits[genome_id][gene] = hit
                elif hit['pident'] > genome_best_hits[genome_id][gene]['pident']:
                    genome_best_hits[genome_id][gene] = hit
    
    print(f"Found {len(genome_best_hits)} genomes with operon hits")
    
    # Extract sequences
    gene_sequences = defaultdict(list)
    processed = 0
    files_not_found = 0
    sequences_extracted = 0
    
    for genome_id, best_hits in genome_best_hits.items():
        processed += 1
        if processed % 100 == 0:
            print(f"  Processed {processed}/{len(genome_best_hits)} genomes")
        
        # Optional: require complete operon (all 7 genes)
        if require_complete and not all(g in best_hits for g in operon_genes):
            continue
        
        # Find genome file
        if source == 'prokka' and prokka_dir:
            genome_file = _find_prokka_fna(prokka_dir, genome_id)
        elif source == 'assemblies' and assemblies_dir:
            genome_file = _find_assembly_fasta(assemblies_dir, genome_id)
        else:
            continue
        
        if not genome_file:
            files_not_found += 1
            if files_not_found <= 10:  # Debug first 10 missing files
                print(f"  Warning: Could not find genome file for {genome_id}")
            continue
        
        # Read genome sequences
        sequences = {}
        records_by_id = _read_genome_sequences(genome_file)
        
        if not records_by_id:
            if processed <= 5:
                print(f"    Warning: No sequences read from {genome_file}")
            continue
        
        for gene, hit in best_hits.items():
            if gene not in operon_genes:
                continue
            
            # Use normalized ID matching for robust contig lookup
            record = None
            hit_vars = _normalize_variants(hit['sseqid'])
            for rid, rec in records_by_id.items():
                if not hit_vars.isdisjoint(_normalize_variants(rid)):
                    record = rec
                    break
            
            if record is None:
                if processed <= 5:  # Debug first few missing contigs
                    print(f"    Warning: Could not find contig {hit['sseqid']} in {genome_id}")
                continue
            
            # Extract sequence based on coordinates
            start = min(hit['sstart'], hit['send']) - 1
            end = max(hit['sstart'], hit['send'])
            
            if hit['sstart'] > hit['send']:
                seq = record.seq[start:end].reverse_complement()
            else:
                seq = record.seq[start:end]
            
            sequences[gene] = seq
        
        # Add sequences to collection
        if sequences:  # Only count if we found any sequences
            sequences_extracted += 1
        for gene in operon_genes:
            if gene in sequences:
                gene_sequences[gene].append({
                    'genome_id': genome_id,
                    'sequence': str(sequences[gene])
                })
    
    # Write sequences to files
    for gene, seqs in gene_sequences.items():
        output_file = os.path.join(output_dir, f"{gene}.fasta")
        with open(output_file, 'w') as f:
            for seq_info in seqs:
                f.write(f">{seq_info['genome_id']}\n")
                f.write(f"{seq_info['sequence']}\n")
        print(f"  Wrote {len(seqs)} sequences for {gene}")
    
    print(f"\nExtraction summary:")
    print(f"  Genomes with BLAST hits: {len(genome_best_hits)}")
    print(f"  Genome files not found: {files_not_found}")
    print(f"  Genomes with extracted sequences: {sequences_extracted}")
    print(f"✅ Extraction complete: {len(gene_sequences)} genes")


def extract_noncoding_sequences(blast_dir: str, assemblies_dir: str,
                               output_dir: str, min_identity: float = 80.0,
                               min_coverage: float = 70.0,
                               source: str = 'assemblies') -> None:
    """
    Extract non-coding sequences (promoter regions) from genomes.
    
    Args:
        blast_dir: Directory with BLAST results for non-coding regions
        assemblies_dir: Directory with genome assemblies
        output_dir: Output directory for sequences
        min_identity: Minimum identity threshold
        min_coverage: Minimum coverage threshold
        source: Source type (assemblies)
    """
    print(f"\n{'='*60}")
    print("Extracting non-coding sequences")
    print(f"{'='*60}")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Process BLAST results for promoter - look for noncoding files
    blast_files = [f for f in os.listdir(blast_dir) if "_noncoding_blast.txt" in f]
    
    if not blast_files:
        print(f"Error: No non-coding BLAST files found in {blast_dir}")
        sys.exit(2)
    
    promoter_sequences = []
    
    for blast_file in blast_files:
        genome_id = _infer_genome_id_from_blast_filename(blast_file)
        file_path = os.path.join(blast_dir, blast_file)
        
        best_hit = None
        with open(file_path, 'r') as f:
            for line in f:
                if not line.strip():
                    continue
                
                hit = parse_blast_hit(line)
                
                if 'promoter' not in hit['qseqid'].lower():
                    continue
                
                if hit['pident'] < min_identity or hit['qcovs'] < min_coverage:
                    continue
                
                if best_hit is None or hit['pident'] > best_hit['pident']:
                    best_hit = hit
        
        if best_hit:
            # Find and extract sequence
            genome_file = _find_assembly_fasta(assemblies_dir, genome_id)
            if genome_file:
                records = _read_genome_sequences(genome_file)
                
                # Try normalized matching
                record = None
                hit_vars = _normalize_variants(best_hit['sseqid'])
                for rid, rec in records.items():
                    if not hit_vars.isdisjoint(_normalize_variants(rid)):
                        record = rec
                        break
                
                if record:
                    start = min(best_hit['sstart'], best_hit['send']) - 1
                    end = max(best_hit['sstart'], best_hit['send'])
                    
                    if best_hit['sstart'] > best_hit['send']:
                        seq = record.seq[start:end].reverse_complement()
                    else:
                        seq = record.seq[start:end]
                    
                    promoter_sequences.append({
                        'genome_id': genome_id,
                        'sequence': str(seq)
                    })
    
    # Write promoter sequences
    if promoter_sequences:
        output_file = os.path.join(output_dir, 'promoter.fasta')
        with open(output_file, 'w') as f:
            for seq_info in promoter_sequences:
                f.write(f">{seq_info['genome_id']}\n")
                f.write(f"{seq_info['sequence']}\n")
        print(f"  Wrote {len(promoter_sequences)} promoter sequences")
    
    print(f"✅ Non-coding extraction complete")


# ============================================================================
# MSA CREATION MODULE
# ============================================================================

def create_msa(sequences_dir: str, output_dir: str, threads: int = 8) -> None:
    """
    Create multiple sequence alignments for DNA sequences.
    
    Args:
        sequences_dir: Directory containing sequence FASTA files
        output_dir: Directory for output alignments
        threads: Number of threads for MAFFT
    """
    print(f"\n{'='*60}")
    print("Creating multiple sequence alignments")
    print(f"{'='*60}")
    
    # Create output directories
    dna_dir = os.path.join(output_dir, 'dna_alignments')
    os.makedirs(dna_dir, exist_ok=True)
    
    # Find sequence files
    seq_files = [f for f in os.listdir(sequences_dir) 
                 if f.endswith('.fasta') or f.endswith('.fa')]
    
    if not seq_files:
        print(f"Warning: No sequence files found in {sequences_dir}")
        return
    
    print(f"Found {len(seq_files)} sequence files to align")
    
    for seq_file in seq_files:
        input_path = os.path.join(sequences_dir, seq_file)
        gene_name = os.path.splitext(seq_file)[0]
        output_path = os.path.join(dna_dir, f"{gene_name}_aligned.fasta")
        
        # Count sequences
        seq_count = sum(1 for line in open(input_path) if line.startswith('>'))
        
        if seq_count < 2:
            print(f"  Skipping {gene_name}: only {seq_count} sequence(s)")
            continue
        
        print(f"  Aligning {gene_name} ({seq_count} sequences)...", end='', flush=True)
        
        # Run MAFFT
        cmd = [
            'mafft',
            '--auto',
            '--thread', str(threads),
            '--quiet',
            input_path
        ]
        
        try:
            with open(output_path, 'w') as out_f:
                result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, 
                                      text=True, check=True)
            print(" ✓")
        except subprocess.CalledProcessError as e:
            print(f" ✗ Error: {e.stderr}")
            continue
    
    print(f"✅ MSA creation complete")


# ============================================================================
# CONSERVATION ANALYSIS MODULE
# ============================================================================

def calculate_shannon_entropy(alignment_file: str) -> Tuple[List[float], Dict]:
    """
    Calculate Shannon entropy-based conservation scores.

    Returns:
        Tuple of (conservation_scores, metadata_dict)
    """
    try:
        sequences = list(SeqIO.parse(alignment_file, 'fasta'))
    except Exception as e:
        print(f"Warning: Could not parse alignment file {alignment_file}: {e}")
        return [], {}

    if not sequences:
        return [], {}

    # Check if all sequences have the same length
    seq_lengths = [len(seq.seq) for seq in sequences]
    if len(set(seq_lengths)) > 1:
        print(f"Warning: Sequences in {alignment_file} have different lengths: {set(seq_lengths)}")
        # Use the minimum length
        alignment_length = min(seq_lengths)
    else:
        alignment_length = len(sequences[0].seq)
    conservation_scores = []
    
    # Variables for fast pairwise identity calculation
    pairwise_identity_sum = 0.0
    pairwise_identity_cols = 0
    
    for pos in range(alignment_length):
        # Count nucleotides at this position
        char_counts = defaultdict(int)
        total_chars = 0
        for seq in sequences:
            if pos < len(seq.seq):
                char = str(seq.seq[pos]).upper()
                if char in 'ACGT':
                    char_counts[char] += 1
                    total_chars += 1
        
        # Calculate Shannon entropy
        if total_chars == 0:
            conservation_scores.append(0.0)
            continue
        
        entropy = 0
        for count in char_counts.values():
            p = count / total_chars
            if p > 0:
                entropy -= p * np.log2(p)
        
        # Convert to conservation score (0-1)
        max_entropy = 2.0  # DNA alphabet - fixed 2-bit normalization
        conservation = 1 - (entropy / max_entropy) if max_entropy > 0 else 1.0
        conservation_scores.append(conservation)
        
        # Fast pairwise identity calculation using the count-based formula
        if total_chars >= 2:
            num_pairs_equal = sum(c*(c-1)//2 for c in char_counts.values())
            den_pairs = total_chars*(total_chars-1)//2
            if den_pairs > 0:
                pairwise_identity_sum += (num_pairs_equal / den_pairs)
                pairwise_identity_cols += 1
    
    # Calculate mean pairwise identity
    mean_pairwise_identity = (pairwise_identity_sum / pairwise_identity_cols * 100) if pairwise_identity_cols > 0 else 0.0
    
    # Calculate metadata
    metadata = {
        'num_sequences': len(sequences),
        'alignment_length': alignment_length,
        'mean_conservation': np.mean(conservation_scores),
        'gap_percentage': calculate_gap_percentage(sequences),
        'pairwise_identity': mean_pairwise_identity
    }
    
    return conservation_scores, metadata


def calculate_detailed_conservation_metrics(alignment_file: str) -> dict:
    """
    Calculate detailed conservation metrics matching core gene analysis format.
    
    Returns dict with conservation metrics including median, std, and position categories.
    """
    from Bio import AlignIO
    from collections import Counter
    
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
            char_counts = Counter([c.upper() for c in column if c != '-'])
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
            'error': str(e),
            'n_sequences': 0
        }


def create_conservation_plots(msa_dir: str, output_dir: str, 
                            title_suffix: str = "") -> None:
    """
    Create enhanced conservation plots with Shannon entropy and sequence logos.
    
    Args:
        msa_dir: Directory containing alignment files
        output_dir: Directory for output plots
        title_suffix: Suffix for plot titles
    """
    print(f"\n{'='*60}")
    print("Creating conservation plots")
    print(f"{'='*60}")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Find alignment files
    align_files = [f for f in os.listdir(msa_dir) 
                   if f.endswith('_aligned.fasta') or f.endswith('.fasta')]
    
    if not align_files:
        print(f"Warning: No alignment files found in {msa_dir}")
        return
    
    for align_file in align_files:
        align_path = os.path.join(msa_dir, align_file)
        gene_name = align_file.replace('_aligned.fasta', '').replace('.fasta', '')
        
        print(f"  Processing {gene_name}...")
        
        # Calculate conservation
        conservation_scores, metadata = calculate_shannon_entropy(align_path)
        
        if not conservation_scores:
            print(f"    Warning: No valid alignment for {gene_name}")
            continue
        
        # Create Shannon entropy plot
        fig, ax = plt.subplots(figsize=(15, 6))
        positions = np.arange(len(conservation_scores))
        
        ax.plot(positions, conservation_scores, color='blue', linewidth=1.5, alpha=0.8)
        ax.fill_between(positions, conservation_scores, alpha=0.3, color='blue')
        
        ax.set_xlabel('Position', fontsize=12)
        ax.set_ylabel('Conservation Score (Shannon Entropy)', fontsize=12)
        ax.set_title(f'{gene_name} Conservation{" - " + title_suffix if title_suffix else ""}',
                    fontsize=14, fontweight='bold')
        ax.set_ylim(0, 1.05)
        ax.grid(True, alpha=0.3)
        
        # Add metadata text
        info_text = (f"Sequences: {metadata['num_sequences']:,}\n"
                    f"Length: {metadata['alignment_length']:,} bp\n"
                    f"Mean conservation: {metadata['mean_conservation']:.3f}")
        ax.text(0.02, 0.98, info_text, transform=ax.transAxes,
                fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        plt.tight_layout()
        output_file = os.path.join(output_dir, f"{gene_name}_shannon_entropy_conservation.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        # Create sequence logo if logomaker available
        if LOGOMAKER_AVAILABLE and metadata['num_sequences'] > 10:
            try:
                create_sequence_logo(align_path, gene_name, output_dir, title_suffix)
            except Exception as e:
                print(f"    Warning: Could not create logo for {gene_name}: {e}")
    
    print(f"✅ Conservation plots complete")


def calculate_gap_fraction(alignment_file: str) -> Tuple[List[float], Dict]:
    """
    Calculate gap fraction for each position in an alignment.
    
    Returns:
        Tuple of (gap_fractions, metadata_dict)
    """
    from Bio import AlignIO
    
    try:
        alignment = AlignIO.read(alignment_file, "fasta")
    except Exception as e:
        print(f"Warning: Could not parse alignment file {alignment_file}: {e}")
        return [], {}
    
    if len(alignment) == 0:
        return [], {}
    
    n_sequences = len(alignment)
    alignment_length = alignment.get_alignment_length()
    gap_fractions = []
    
    for pos in range(alignment_length):
        column = alignment[:, pos]
        gap_count = column.count('-')
        gap_fraction = gap_count / n_sequences if n_sequences > 0 else 0
        gap_fractions.append(gap_fraction)
    
    # Calculate metadata
    metadata = {
        'num_sequences': n_sequences,
        'alignment_length': alignment_length,
        'mean_gap_fraction': np.mean(gap_fractions),
        'max_gap_fraction': np.max(gap_fractions) if gap_fractions else 0,
        'positions_with_gaps': sum(1 for g in gap_fractions if g > 0),
        'positions_with_gaps_pct': (sum(1 for g in gap_fractions if g > 0) / len(gap_fractions) * 100) if gap_fractions else 0
    }
    
    return gap_fractions, metadata


def create_gap_plots(msa_dir: str, output_dir: str, title_suffix: str = "") -> None:
    """
    Create gap fraction plots for alignments.
    
    Args:
        msa_dir: Directory containing alignment files
        output_dir: Directory for output plots
        title_suffix: Suffix for plot titles
    """
    print(f"\n{'='*60}")
    print("Creating gap fraction plots")
    print(f"{'='*60}")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Find alignment files
    align_files = sorted([f for f in os.listdir(msa_dir) 
                         if f.endswith('_aligned.fasta') or f.endswith('.fasta')])
    
    if not align_files:
        print(f"Warning: No alignment files found in {msa_dir}")
        return
    
    for align_file in align_files:
        align_path = os.path.join(msa_dir, align_file)
        gene_name = align_file.replace('_aligned.fasta', '').replace('.fasta', '')
        
        print(f"  Processing {gene_name}...")
        
        # Calculate gap fractions
        gap_fractions, metadata = calculate_gap_fraction(align_path)
        
        if not gap_fractions:
            print(f"    Warning: No valid alignment for {gene_name}")
            continue
        
        # Create single gap fraction plot
        fig, ax = plt.subplots(figsize=(15, 6))
        positions = np.arange(len(gap_fractions))
        
        # Gap fraction plot without fill overlay
        ax.plot(positions, gap_fractions, color='red', linewidth=1.5, alpha=0.8)
        ax.set_xlabel('Position', fontsize=12)
        ax.set_ylabel('Gap Fraction', fontsize=12)
        ax.set_title(f'{gene_name} Gap Distribution{" - " + title_suffix if title_suffix else ""}',
                    fontsize=14, fontweight='bold')
        ax.set_ylim(0, 1.05)
        ax.grid(True, alpha=0.3)
        
        # Add horizontal lines for reference
        ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, label='50% gaps')
        ax.axhline(y=0.9, color='gray', linestyle=':', alpha=0.5, label='90% gaps')
        
        # Add metadata text
        info_text = (f"Sequences: {metadata['num_sequences']:,}\n"
                    f"Length: {metadata['alignment_length']:,} bp\n"
                    f"Mean gap fraction: {metadata['mean_gap_fraction']:.3f}\n"
                    f"Max gap fraction: {metadata['max_gap_fraction']:.3f}\n"
                    f"Positions with gaps: {metadata['positions_with_gaps_pct']:.1f}%")
        ax.text(0.02, 0.98, info_text, transform=ax.transAxes,
               fontsize=10, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        ax.legend(loc='upper right', fontsize=10)
        
        plt.tight_layout()
        output_file = os.path.join(output_dir, f"{gene_name}_gap_fraction_analysis.png")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"    Saved gap plot to {output_file}")
    
    print(f"✅ Gap fraction plots complete")


def create_sequence_logo(alignment_file: str, gene_name: str,
                        output_dir: str, title_suffix: str = "") -> None:
    """Create sequence logo visualization."""
    if not LOGOMAKER_AVAILABLE:
        print(f"  Skipping sequence logo for {gene_name}: logomaker not available")
        return

    sequences = list(SeqIO.parse(alignment_file, 'fasta'))
    
    if len(sequences) < 10:
        return
    
    # Sample positions for visualization (first 100 bp or full length)
    max_pos = min(100, len(sequences[0].seq))
    
    # Create position frequency matrix
    matrix_data = []
    for pos in range(max_pos):
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        total = 0
        for seq in sequences:
            if pos < len(seq.seq):
                nuc = str(seq.seq[pos]).upper()
                if nuc in counts:
                    counts[nuc] += 1
                    total += 1
        
        # Convert to frequencies
        if total > 0:
            matrix_data.append({n: c/total for n, c in counts.items()})
        else:
            matrix_data.append({'A': 0, 'C': 0, 'G': 0, 'T': 0})
    
    # Create DataFrame for logomaker
    df = pd.DataFrame(matrix_data)
    
    # Create logo
    fig, ax = plt.subplots(figsize=(15, 4))
    logo = logomaker.Logo(df, ax=ax, color_scheme='classic')
    
    ax.set_xlabel('Position', fontsize=12)
    ax.set_ylabel('Frequency', fontsize=12)
    ax.set_title(f'{gene_name} Sequence Logo (first {max_pos} bp){" - " + title_suffix if title_suffix else ""}',
                fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, f"{gene_name}_nucleotide_frequency_logo.png")
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()


def create_msa_summary(msa_dir: str, output_file: str, 
                      strategy_name: str = "Primary Analysis") -> None:
    """
    Create summary CSV with conservation metrics for all alignments.
    
    Args:
        msa_dir: Directory containing MSA subdirectories
        output_file: Output CSV file path
        strategy_name: Name of the analysis strategy
    """
    print(f"\n{'='*60}")
    print("Creating MSA summary metrics")
    print(f"{'='*60}")
    
    # Find DNA alignments
    dna_dir = os.path.join(msa_dir, 'dna_alignments')
    
    if not os.path.exists(dna_dir):
        print(f"Warning: DNA alignment directory not found: {dna_dir}")
        return
    
    results = []
    
    align_files = [f for f in os.listdir(dna_dir) 
                   if f.endswith('_aligned.fasta') or f.endswith('.fasta')]
    
    for align_file in align_files:
        align_path = os.path.join(dna_dir, align_file)
        gene_name = align_file.replace('_aligned.fasta', '').replace('.fasta', '')
        
        print(f"  Processing {gene_name}...")
        
        conservation_scores, metadata = calculate_shannon_entropy(align_path)
        
        if conservation_scores:
            results.append({
                'gene': gene_name,
                'num_sequences': metadata['num_sequences'],
                'alignment_length': metadata['alignment_length'],
                'conservation_score': metadata['mean_conservation'],
                'gap_percentage': metadata['gap_percentage'],
                'pairwise_identity': metadata['pairwise_identity']
            })
    
    # Save to CSV
    if results:
        df = pd.DataFrame(results)
        df = df.sort_values('conservation_score', ascending=False)
        df.to_csv(output_file, index=False)
        print(f"  Saved metrics for {len(results)} genes to {output_file}")
    
    print(f"✅ MSA summary complete")


def create_detailed_conservation_summary(msa_dir: str, output_file: str) -> None:
    """
    Create detailed conservation metrics CSV matching core gene analysis format.
    
    Args:
        msa_dir: Directory containing MSA subdirectories
        output_file: Output CSV file path
    """
    print(f"\n{'='*60}")
    print("Creating detailed conservation metrics")
    print(f"{'='*60}")
    
    # Find DNA alignments
    dna_dir = os.path.join(msa_dir, 'dna_alignments')
    
    if not os.path.exists(dna_dir):
        print(f"Warning: DNA alignment directory not found: {dna_dir}")
        return
    
    results = []
    
    align_files = sorted([f for f in os.listdir(dna_dir) 
                         if f.endswith('_aligned.fasta') or f.endswith('.fasta')])
    
    for i, align_file in enumerate(align_files, 1):
        align_path = os.path.join(dna_dir, align_file)
        gene_name = align_file.replace('_aligned.fasta', '').replace('.fasta', '')
        
        print(f"  Processing {gene_name} ({i}/{len(align_files)})...")
        
        metrics = calculate_detailed_conservation_metrics(align_path)
        
        if 'error' not in metrics:
            results.append(metrics)
        else:
            print(f"    Warning: {metrics.get('error', 'Unknown error')}")
    
    # Save to CSV
    if results:
        df = pd.DataFrame(results)
        df = df.sort_values('mean_conservation', ascending=False)
        df.to_csv(output_file, index=False)
        
        print(f"\n✅ Saved detailed metrics for {len(results)} genes to {output_file}")
        
        # Print summary
        print("\nSummary Statistics:")
        print(f"  Mean conservation: {df['mean_conservation'].mean():.4f} ± {df['mean_conservation'].std():.4f}")
        print(f"  Highly conserved positions: {df['highly_conserved_pct'].mean():.1f}%")
        print(f"  Moderately conserved: {df['moderately_conserved_pct'].mean():.1f}%")
        print(f"  Variable positions: {df['variable_pct'].mean():.1f}%")
        print(f"  Positions with gaps: {df['positions_with_gaps_pct'].mean():.1f}%")
    else:
        print("No successful alignments processed")


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def _find_prokka_fna(prokka_dir: str, genome_id: str) -> Optional[str]:
    """Find Prokka FNA file for a genome."""
    # try exact subdirs
    cand_dirs = [
        os.path.join(prokka_dir, f"{genome_id}.result"),
        os.path.join(prokka_dir, genome_id),
    ]
    # also scan for unique "*.result" dirs prefixed by genome_id
    try:
        for entry in os.listdir(prokka_dir):
            full = os.path.join(prokka_dir, entry)
            if os.path.isdir(full) and entry.startswith(genome_id) and entry.endswith(".result"):
                if full not in cand_dirs:
                    cand_dirs.insert(0, full)
    except FileNotFoundError:
        pass
    for d in cand_dirs:
        if not os.path.isdir(d):
            continue
        base = os.path.basename(d)
        for cand in (os.path.join(d, f"{base}.fna"), os.path.join(d, f"{genome_id}.fna")):
            if os.path.exists(cand):
                return cand
        fna_files = [os.path.join(d, f) for f in os.listdir(d) if f.endswith(".fna")]
        if len(fna_files) == 1:
            return fna_files[0]
        for f in fna_files:
            if genome_id in os.path.basename(f):
                return f
        if fna_files:
            return fna_files[0]
    top_level = os.path.join(prokka_dir, f"{genome_id}.fna")
    return top_level if os.path.exists(top_level) else None


def _find_assembly_fasta(assemblies_dir: str, genome_id: str) -> Optional[str]:
    """Find assembly FASTA file for a genome."""
    # Direct match with standard suffixes
    for suf in (".fasta.gz", ".fasta", ".fa.gz", ".fa"):
        p = os.path.join(assemblies_dir, genome_id + suf)
        if os.path.exists(p):
            return p
    
    # Try adding .result or .genome before fasta extension
    # For files like ENT_AA0002AA_AS.result.fasta.gz or ENT_CA1913AA_AS.genome.fasta.gz
    for middle in (".result", ".genome"):
        for suf in (".fasta.gz", ".fasta", ".fa.gz", ".fa"):
            p = os.path.join(assemblies_dir, genome_id + middle + suf)
            if os.path.exists(p):
                return p
    
    return None


def _read_genome_sequences(genome_file: str) -> Dict:
    """Read sequences from genome file."""
    sequences = {}

    try:
        if genome_file.endswith('.gz'):
            handle = gzip.open(genome_file, 'rt')
        else:
            handle = open(genome_file, 'r')
    except (FileNotFoundError, IOError) as e:
        print(f"Warning: Could not open genome file {genome_file}: {e}")
        return sequences

    try:
        for record in SeqIO.parse(handle, "fasta"):
            sequences[record.id] = record
    except Exception as e:
        print(f"Warning: Error parsing sequences from {genome_file}: {e}")
    finally:
        handle.close()

    return sequences


def calculate_gap_percentage(sequences: List) -> float:
    """Calculate percentage of gaps in alignment."""
    if not sequences:
        return 0.0

    alignment_length = len(sequences[0].seq)
    column_gap_percentages = []

    for pos in range(alignment_length):
        gap_count = sum(1 for seq in sequences if pos < len(seq.seq) and str(seq.seq[pos]) == '-')
        gap_percentage = (gap_count / len(sequences)) * 100
        column_gap_percentages.append(gap_percentage)

    return np.mean(column_gap_percentages) if column_gap_percentages else 0.0


def calculate_pairwise_identity(sequences: List) -> float:
    """Calculate mean pairwise identity using fast count-based method."""
    if len(sequences) < 2:
        return 100.0
    
    # Get alignment length (assuming all sequences are same length)
    alignment_length = min(len(seq.seq) for seq in sequences)
    
    pairwise_identity_sum = 0.0
    pairwise_identity_cols = 0
    
    for pos in range(alignment_length):
        # Count characters at this position (excluding gaps)
        char_counts = defaultdict(int)
        total_chars = 0
        
        for seq in sequences:
            if pos < len(seq.seq):
                char = str(seq.seq[pos]).upper()
                if char != '-' and char in 'ACGT':
                    char_counts[char] += 1
                    total_chars += 1
        
        # Fast pairwise identity calculation using the count-based formula
        if total_chars >= 2:
            num_pairs_equal = sum(c*(c-1)//2 for c in char_counts.values())
            den_pairs = total_chars*(total_chars-1)//2
            if den_pairs > 0:
                pairwise_identity_sum += (num_pairs_equal / den_pairs)
                pairwise_identity_cols += 1
    
    return (pairwise_identity_sum / pairwise_identity_cols * 100) if pairwise_identity_cols > 0 else 0.0


# ============================================================================
# MAIN COMMAND-LINE INTERFACE
# ============================================================================

def main():
    """Main entry point for the consolidated pipeline."""
    parser = argparse.ArgumentParser(
        description='Consolidated operon extraction and analysis pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Commands:
  extract-sequences    Extract operon sequences from assemblies
  extract-noncoding   Extract non-coding sequences (promoters)
  create-msa          Create multiple sequence alignments
  create-plots        Create conservation plots
  create-gap-plots    Create gap fraction plots
  create-summary      Create conservation metrics summary
  create-detailed-summary  Create detailed conservation metrics (matching core gene format)
  
Examples:
  %(prog)s extract-sequences --blast-dir ../03_blast_search/output/blast_results \\
                           --assemblies-dir ../../Efs_assemblies \\
                           --output-dir output/sequences
  
  %(prog)s create-msa --sequences-dir output/sequences \\
                     --output-dir output/msa
  
  %(prog)s create-plots --msa-dir output/msa/dna_alignments \\
                       --output-dir output/plots
""")
    
    subparsers = parser.add_subparsers(dest='command', help='Command to run')
    
    # Extract sequences command
    extract_parser = subparsers.add_parser('extract-sequences',
                                          help='Extract operon sequences')
    extract_parser.add_argument('--blast-dir', required=True,
                               help='Directory with BLAST results')
    extract_parser.add_argument('--output-dir', required=True,
                               help='Output directory for sequences')
    extract_parser.add_argument('--prokka-dir',
                               help='Prokka results directory (for source=prokka)')
    extract_parser.add_argument('--assemblies-dir',
                               help='Assemblies directory (for source=assemblies)')
    extract_parser.add_argument('--min-identity', type=float, default=90.0,
                               help='Minimum identity threshold (default: 90)')
    extract_parser.add_argument('--min-coverage', type=float, default=80.0,
                               help='Minimum coverage threshold (default: 80)')
    extract_parser.add_argument('--source', choices=['prokka', 'assemblies'],
                               default='assemblies',
                               help='Source for sequences (default: assemblies)')
    extract_parser.add_argument('--require-complete', action='store_true',
                               help='Only extract from genomes with all 7 operon genes')
    
    # Extract non-coding command
    noncoding_parser = subparsers.add_parser('extract-noncoding',
                                            help='Extract non-coding sequences')
    noncoding_parser.add_argument('--blast-dir', required=True,
                                 help='Directory with BLAST results')
    noncoding_parser.add_argument('--assemblies-dir', required=True,
                                 help='Assemblies directory')
    noncoding_parser.add_argument('--output-dir', required=True,
                                 help='Output directory')
    noncoding_parser.add_argument('--min-identity', type=float, default=80.0,
                                 help='Minimum identity (default: 80)')
    noncoding_parser.add_argument('--min-coverage', type=float, default=70.0,
                                 help='Minimum coverage (default: 70)')
    noncoding_parser.add_argument('--source', default='assemblies',
                                 help='Source type (default: assemblies)')
    
    # Create MSA command
    msa_parser = subparsers.add_parser('create-msa',
                                       help='Create multiple sequence alignments')
    msa_parser.add_argument('--sequences-dir', required=True,
                           help='Directory with sequence files')
    msa_parser.add_argument('--output-dir', required=True,
                           help='Output directory for alignments')
    msa_parser.add_argument('--threads', type=int, default=8,
                           help='Number of threads for MAFFT (default: 8)')
    
    # Create plots command
    plots_parser = subparsers.add_parser('create-plots',
                                        help='Create conservation plots')
    plots_parser.add_argument('--msa-dir', required=True,
                             help='Directory with alignment files')
    plots_parser.add_argument('--output-dir', required=True,
                             help='Output directory for plots')
    plots_parser.add_argument('--title-suffix', default='',
                             help='Suffix for plot titles')
    
    # Create gap plots command
    gap_plots_parser = subparsers.add_parser('create-gap-plots',
                                            help='Create gap fraction plots')
    gap_plots_parser.add_argument('--msa-dir', required=True,
                                 help='Directory with alignment files')
    gap_plots_parser.add_argument('--output-dir', required=True,
                                 help='Output directory for plots')
    gap_plots_parser.add_argument('--title-suffix', default='',
                                 help='Suffix for plot titles')
    
    # Create summary command
    summary_parser = subparsers.add_parser('create-summary',
                                          help='Create conservation summary')
    summary_parser.add_argument('--msa-dir', required=True,
                               help='Directory with MSA subdirectories')
    summary_parser.add_argument('--output-file', required=True,
                               help='Output CSV file')
    summary_parser.add_argument('--strategy-name', default='Primary Analysis',
                               help='Analysis strategy name')
    
    # Create detailed summary command
    detailed_parser = subparsers.add_parser('create-detailed-summary',
                                           help='Create detailed conservation metrics matching core gene format')
    detailed_parser.add_argument('--msa-dir', required=True,
                                help='Directory with MSA subdirectories')
    detailed_parser.add_argument('--output-file', required=True,
                                help='Output CSV file')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        sys.exit(1)
    
    # Execute command
    if args.command == 'extract-sequences':
        extract_operon_sequences(
            blast_dir=args.blast_dir,
            output_dir=args.output_dir,
            prokka_dir=args.prokka_dir,
            assemblies_dir=args.assemblies_dir,
            min_identity=args.min_identity,
            min_coverage=args.min_coverage,
            source=args.source,
            require_complete=getattr(args, 'require_complete', False)
        )
    
    elif args.command == 'extract-noncoding':
        # Use the main blast_results directory for noncoding files
        blast_dir = args.blast_dir
        if "blast_results_noncoding" in blast_dir:
            blast_dir = blast_dir.replace("blast_results_noncoding", "blast_results")
        
        extract_noncoding_sequences(
            blast_dir=blast_dir,
            assemblies_dir=args.assemblies_dir,
            output_dir=args.output_dir,
            min_identity=args.min_identity,
            min_coverage=args.min_coverage,
            source=args.source
        )
    
    elif args.command == 'create-msa':
        create_msa(
            sequences_dir=args.sequences_dir,
            output_dir=args.output_dir,
            threads=args.threads
        )
    
    elif args.command == 'create-plots':
        create_conservation_plots(
            msa_dir=args.msa_dir,
            output_dir=args.output_dir,
            title_suffix=args.title_suffix
        )
    
    elif args.command == 'create-gap-plots':
        create_gap_plots(
            msa_dir=args.msa_dir,
            output_dir=args.output_dir,
            title_suffix=args.title_suffix
        )
    
    elif args.command == 'create-summary':
        create_msa_summary(
            msa_dir=args.msa_dir,
            output_file=args.output_file,
            strategy_name=args.strategy_name
        )
    
    elif args.command == 'create-detailed-summary':
        create_detailed_conservation_summary(
            msa_dir=args.msa_dir,
            output_file=args.output_file
        )


if __name__ == '__main__':
    main()