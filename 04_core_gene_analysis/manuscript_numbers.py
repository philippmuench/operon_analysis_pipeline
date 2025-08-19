#!/usr/bin/env python3
"""
Generate manuscript statistics for core gene analysis methods section.
"""

import os
import pandas as pd
import glob
from Bio import SeqIO

def analyze_core_gene_identification():
    """Analyze core gene identification results."""
    stats = {}
    
    # Check for core genes file
    core_genes_file = "output/core_genes_95pct.txt"
    if os.path.exists(core_genes_file):
        with open(core_genes_file, 'r') as f:
            core_genes = [line.strip() for line in f if line.strip()]
        stats['core_genes_count'] = len(core_genes)
    else:
        stats['core_genes_count'] = 0
    
    # Check for prevalence statistics
    prevalence_file = "output/gene_prevalence_stats.csv"
    if os.path.exists(prevalence_file):
        try:
            df = pd.read_csv(prevalence_file)
            stats['total_unique_genes'] = len(df)
            stats['prevalence_distribution'] = {}
            
            # Calculate prevalence thresholds
            for threshold in [1.0, 0.99, 0.95, 0.90, 0.80, 0.50]:
                count = (df['prevalence'] >= threshold).sum()
                stats['prevalence_distribution'][f'{threshold*100:.0f}%'] = count
            
            # Additional statistics
            stats['mean_prevalence'] = df['prevalence'].mean()
            stats['genes_in_all_genomes'] = (df['prevalence'] == 1.0).sum()
            stats['rare_genes_10pct'] = (df['prevalence'] < 0.1).sum()
            
        except Exception as e:
            stats['prevalence_error'] = str(e)
    
    return stats

def analyze_sequence_extraction():
    """Analyze sequence extraction results."""
    stats = {}
    
    sequences_dir = "output/core_gene_sequences"
    if os.path.exists(sequences_dir):
        fasta_files = glob.glob(os.path.join(sequences_dir, "*.fasta"))
        stats['genes_extracted'] = len(fasta_files)
        
        # Analyze a sample of sequences
        total_sequences = 0
        sequence_lengths = []
        
        for fasta_file in fasta_files[:10]:  # Sample first 10 files
            try:
                with open(fasta_file, 'r') as f:
                    count = 0
                    for record in SeqIO.parse(f, 'fasta'):
                        count += 1
                        sequence_lengths.append(len(record.seq))
                    total_sequences += count
            except:
                continue
        
        if sequence_lengths:
            stats['avg_sequences_per_gene'] = total_sequences / min(10, len(fasta_files))
            stats['avg_sequence_length'] = sum(sequence_lengths) / len(sequence_lengths)
            stats['min_sequence_length'] = min(sequence_lengths)
            stats['max_sequence_length'] = max(sequence_lengths)
    else:
        stats['genes_extracted'] = 0
    
    return stats

def analyze_msa_creation():
    """Analyze MSA creation results."""
    stats = {}
    
    alignments_dir = "output/core_gene_alignments"
    if os.path.exists(alignments_dir):
        alignment_files = glob.glob(os.path.join(alignments_dir, "*_aligned.fasta"))
        stats['alignments_created'] = len(alignment_files)
        
        # Analyze alignment properties
        alignment_lengths = []
        sequence_counts = []
        
        for align_file in alignment_files[:10]:  # Sample first 10
            try:
                with open(align_file, 'r') as f:
                    records = list(SeqIO.parse(f, 'fasta'))
                    if records:
                        alignment_lengths.append(len(records[0].seq))
                        sequence_counts.append(len(records))
            except:
                continue
        
        if alignment_lengths:
            stats['avg_alignment_length'] = sum(alignment_lengths) / len(alignment_lengths)
            stats['avg_sequences_per_alignment'] = sum(sequence_counts) / len(sequence_counts)
    else:
        stats['alignments_created'] = 0
    
    return stats

def analyze_conservation_metrics():
    """Analyze conservation analysis results."""
    stats = {}
    
    conservation_file = "output/core_gene_conservation_metrics.csv"
    if os.path.exists(conservation_file):
        try:
            df = pd.read_csv(conservation_file)
            stats['genes_analyzed'] = len(df)
            stats['mean_conservation_score'] = df['mean_conservation'].mean()
            stats['std_conservation_score'] = df['mean_conservation'].std()
            
            # Conservation categories
            if 'conservation_category' in df.columns:
                category_counts = df['conservation_category'].value_counts()
                stats['conservation_categories'] = category_counts.to_dict()
            
            # Identity metrics
            if 'mean_pairwise_identity' in df.columns:
                stats['mean_pairwise_identity'] = df['mean_pairwise_identity'].mean()
            
            # Gap statistics
            if 'gap_percentage' in df.columns:
                stats['mean_gap_percentage'] = df['gap_percentage'].mean()
                
        except Exception as e:
            stats['conservation_error'] = str(e)
    
    return stats

def count_input_genomes():
    """Count input genomes from Prokka results."""
    prokka_dir = "../01_prokka_annotation/output/prokka_results"
    if os.path.exists(prokka_dir):
        genome_dirs = glob.glob(os.path.join(prokka_dir, "*/"))
        return len(genome_dirs)
    return 0

def generate_manuscript_stats():
    """Generate all statistics for manuscript."""
    
    print("Core Gene Analysis Statistics for Manuscript")
    print("=" * 60)
    
    # Input genomes
    input_genomes = count_input_genomes()
    print(f"\n1. Input Data:")
    print(f"   Genomes analyzed: {input_genomes:,}")
    
    # Core gene identification
    print(f"\n2. Core Gene Identification:")
    core_stats = analyze_core_gene_identification()
    print(f"   Total unique genes: {core_stats.get('total_unique_genes', 0):,}")
    print(f"   Core genes (≥95%): {core_stats.get('core_genes_count', 0):,}")
    print(f"   Genes in 100% genomes: {core_stats.get('genes_in_all_genomes', 0):,}")
    
    if 'prevalence_distribution' in core_stats:
        print(f"   Prevalence distribution:")
        for threshold, count in core_stats['prevalence_distribution'].items():
            print(f"     ≥{threshold}: {count:,} genes")
    
    # Sequence extraction
    print(f"\n3. Sequence Extraction:")
    seq_stats = analyze_sequence_extraction()
    print(f"   Genes with sequences extracted: {seq_stats.get('genes_extracted', 0):,}")
    if 'avg_sequences_per_gene' in seq_stats:
        print(f"   Average sequences per gene: {seq_stats['avg_sequences_per_gene']:.0f}")
        print(f"   Average sequence length: {seq_stats['avg_sequence_length']:.0f} bp")
        print(f"   Length range: {seq_stats['min_sequence_length']}-{seq_stats['max_sequence_length']} bp")
    
    # MSA creation
    print(f"\n4. Multiple Sequence Alignments:")
    msa_stats = analyze_msa_creation()
    print(f"   Alignments created: {msa_stats.get('alignments_created', 0):,}")
    if 'avg_alignment_length' in msa_stats:
        print(f"   Average alignment length: {msa_stats['avg_alignment_length']:.0f} positions")
        print(f"   Average sequences per alignment: {msa_stats['avg_sequences_per_alignment']:.0f}")
    
    # Conservation analysis
    print(f"\n5. Conservation Analysis:")
    cons_stats = analyze_conservation_metrics()
    if 'genes_analyzed' in cons_stats:
        print(f"   Genes analyzed: {cons_stats['genes_analyzed']:,}")
        print(f"   Mean conservation score: {cons_stats['mean_conservation_score']:.3f}")
        if 'conservation_categories' in cons_stats:
            print(f"   Conservation categories:")
            for category, count in cons_stats['conservation_categories'].items():
                print(f"     {category}: {count:,} genes")
        
        if 'mean_pairwise_identity' in cons_stats:
            print(f"   Mean pairwise identity: {cons_stats['mean_pairwise_identity']:.3f}")
        if 'mean_gap_percentage' in cons_stats:
            print(f"   Mean gap percentage: {cons_stats['mean_gap_percentage']:.1f}%")
    else:
        print(f"   Conservation analysis not completed")
    
    # Summary for manuscript
    print("\n" + "=" * 60)
    print("MANUSCRIPT NUMBERS SUMMARY:")
    print("=" * 60)
    print(f"Input genomes: {input_genomes:,}")
    print(f"Unique genes identified: {core_stats.get('total_unique_genes', 0):,}")
    print(f"Core genes (≥95% prevalence): {core_stats.get('core_genes_count', 0):,}")
    print(f"Genes in 100% of genomes: {core_stats.get('genes_in_all_genomes', 0):,}")
    print(f"Core gene sequences extracted: {seq_stats.get('genes_extracted', 0):,}")
    print(f"Multiple sequence alignments: {msa_stats.get('alignments_created', 0):,}")
    if 'mean_conservation_score' in cons_stats:
        print(f"Mean conservation score: {cons_stats['mean_conservation_score']:.3f}")
    print(f"Analysis threshold: ≥95% prevalence")
    print(f"Software: MAFFT for alignments, Shannon entropy for conservation")
    
    return {
        'input_genomes': input_genomes,
        'core_stats': core_stats,
        'seq_stats': seq_stats,
        'msa_stats': msa_stats,
        'cons_stats': cons_stats
    }

if __name__ == "__main__":
    # Change to script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    
    stats = generate_manuscript_stats()
