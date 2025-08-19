#!/usr/bin/env python3
"""
Generate manuscript statistics for BLAST search methods section.
"""

import os
import pandas as pd
import glob
from collections import defaultdict

def count_blast_files_and_hits():
    """Count BLAST files and hits across all search modes."""
    stats = {}
    
    # Count files and hits for each search mode
    search_modes = [
        ("coding_protein", "output/blast_results", "*_genes_blast.txt"),
        ("coding_nt", "output/blast_results_nt", "*_genes_blast.txt"),
        ("prokka_variants", "output/blast_results_prokka_variants", "*_blast.txt"),
        ("noncoding", "output/blast_results", "*_noncoding_blast.txt")
    ]
    
    for mode_name, directory, pattern in search_modes:
        if not os.path.exists(directory):
            stats[mode_name] = {"files": 0, "total_hits": 0, "genomes": 0}
            continue
            
        files = glob.glob(os.path.join(directory, pattern))
        total_hits = 0
        non_empty_files = 0
        
        for file_path in files:
            if os.path.getsize(file_path) > 0:
                non_empty_files += 1
                try:
                    with open(file_path, 'r') as f:
                        hits = sum(1 for line in f if line.strip())
                    total_hits += hits
                except:
                    continue
        
        stats[mode_name] = {
            "files": len(files),
            "non_empty_files": non_empty_files,
            "total_hits": total_hits,
            "genomes": non_empty_files
        }
    
    return stats

def analyze_hit_quality():
    """Analyze hit quality distribution."""
    quality_stats = {}
    
    # Analyze main coding protein results
    blast_dir = "output/blast_results"
    if os.path.exists(blast_dir):
        gene_files = glob.glob(os.path.join(blast_dir, '*_genes_blast.txt'))
        
        all_identities = []
        all_coverages = []
        high_quality_hits = 0
        total_hits = 0
        
        for file_path in gene_files:
            if os.path.getsize(file_path) == 0:
                continue
                
            try:
                df = pd.read_csv(file_path, sep='\t', header=None,
                               names=['qseqid', 'sseqid', 'pident', 'length', 
                                     'mismatch', 'gapopen', 'qstart', 'qend',
                                     'sstart', 'send', 'evalue', 'bitscore', 'qcovs'])
                
                all_identities.extend(df['pident'].tolist())
                all_coverages.extend(df['qcovs'].tolist())
                
                # Count high-quality hits (≥90% identity, ≥80% coverage)
                high_quality = df[(df['pident'] >= 90) & (df['qcovs'] >= 80)]
                high_quality_hits += len(high_quality)
                total_hits += len(df)
                
            except Exception as e:
                continue
        
        quality_stats = {
            "total_hits": total_hits,
            "high_quality_hits": high_quality_hits,
            "high_quality_percent": (high_quality_hits / total_hits * 100) if total_hits > 0 else 0,
            "mean_identity": sum(all_identities) / len(all_identities) if all_identities else 0,
            "mean_coverage": sum(all_coverages) / len(all_coverages) if all_coverages else 0
        }
    
    return quality_stats

def count_reference_sequences():
    """Count reference sequences used as queries."""
    ref_counts = {}
    
    ref_files = [
        ("operon_genes_protein", "../02_reference_operon_extraction/output/operon_genes_protein.fasta"),
        ("operon_genes_nt", "../02_reference_operon_extraction/output/operon_genes_nt.fasta"),
        ("operon_noncoding_nt", "../02_reference_operon_extraction/output/operon_noncoding_nt.fasta")
    ]
    
    for name, file_path in ref_files:
        if os.path.exists(file_path):
            try:
                with open(file_path, 'r') as f:
                    count = sum(1 for line in f if line.startswith('>'))
                ref_counts[name] = count
            except:
                ref_counts[name] = 0
        else:
            ref_counts[name] = 0
    
    return ref_counts

def analyze_operon_completeness():
    """Analyze operon completeness across genomes."""
    completeness_stats = {}
    
    # Check if summary file exists
    summary_file = "output/operon_simple_summary.csv"
    if os.path.exists(summary_file):
        try:
            df = pd.read_csv(summary_file)
            
            completeness_stats = {
                "total_genomes": len(df),
                "complete_operons": df['is_complete'].sum(),
                "complete_percent": (df['is_complete'].sum() / len(df) * 100),
                "has_promoter": df['has_promoter'].sum() if 'has_promoter' in df.columns else 0,
                "mean_completeness": df['completeness_score'].mean() if 'completeness_score' in df.columns else 0
            }
        except:
            completeness_stats = {"error": "Could not read summary file"}
    
    return completeness_stats

def generate_manuscript_stats():
    """Generate all statistics for manuscript."""
    
    print("BLAST Search Statistics for Manuscript")
    print("=" * 50)
    
    # Reference sequences
    print("\n1. Reference Query Sequences:")
    ref_counts = count_reference_sequences()
    for name, count in ref_counts.items():
        print(f"   {name}: {count} sequences")
    
    # BLAST search results
    print("\n2. BLAST Search Results:")
    blast_stats = count_blast_files_and_hits()
    for mode, stats in blast_stats.items():
        print(f"   {mode}:")
        print(f"     - Genomes searched: {stats['genomes']}")
        print(f"     - Total hits: {stats['total_hits']:,}")
    
    # Hit quality
    print("\n3. Hit Quality Analysis:")
    quality_stats = analyze_hit_quality()
    if quality_stats:
        print(f"   Total hits analyzed: {quality_stats['total_hits']:,}")
        print(f"   High-quality hits (≥90% identity, ≥80% coverage): {quality_stats['high_quality_hits']:,} ({quality_stats['high_quality_percent']:.1f}%)")
        print(f"   Mean sequence identity: {quality_stats['mean_identity']:.1f}%")
        print(f"   Mean query coverage: {quality_stats['mean_coverage']:.1f}%")
    
    # Operon completeness
    print("\n4. Operon Completeness:")
    completeness_stats = analyze_operon_completeness()
    if "total_genomes" in completeness_stats:
        print(f"   Total genomes analyzed: {completeness_stats['total_genomes']:,}")
        print(f"   Complete operons: {completeness_stats['complete_operons']:,} ({completeness_stats['complete_percent']:.1f}%)")
        if completeness_stats['has_promoter'] > 0:
            print(f"   Genomes with promoter: {completeness_stats['has_promoter']:,}")
    
    # Summary for manuscript
    print("\n" + "=" * 50)
    print("MANUSCRIPT NUMBERS SUMMARY:")
    print("=" * 50)
    
    total_protein_queries = ref_counts.get('operon_genes_protein', 0)
    total_nt_queries = ref_counts.get('operon_genes_nt', 0) + ref_counts.get('operon_noncoding_nt', 0)
    total_genomes = blast_stats.get('coding_protein', {}).get('genomes', 0)
    total_hits = sum(stats.get('total_hits', 0) for stats in blast_stats.values())
    
    print(f"Query sequences: {total_protein_queries} protein, {total_nt_queries} nucleotide")
    print(f"Genomes searched: {total_genomes:,}")
    print(f"Total BLAST hits: {total_hits:,}")
    if quality_stats:
        print(f"High-quality hits: {quality_stats['high_quality_hits']:,} ({quality_stats['high_quality_percent']:.1f}%)")
    
    return {
        'reference_counts': ref_counts,
        'blast_stats': blast_stats,
        'quality_stats': quality_stats,
        'completeness_stats': completeness_stats
    }

if __name__ == "__main__":
    # Change to script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    
    stats = generate_manuscript_stats()
