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

def generate_manuscript_stats(output_file=None):
    """Generate all statistics for manuscript.
    
    Args:
        output_file (str): Optional path to save output to file. If None, prints to console.
    """
    
    # Collect all output lines
    output_lines = []
    
    def add_line(text=""):
        """Add a line to both console and file output."""
        print(text)
        output_lines.append(text)
    
    add_line("BLAST Search Statistics for Manuscript")
    add_line("=" * 60)
    
    # Reference sequences
    add_line("\n1. Reference Query Sequences:")
    ref_counts = count_reference_sequences()
    total_protein = ref_counts.get('operon_genes_protein', 0)
    total_nt_genes = ref_counts.get('operon_genes_nt', 0)
    total_noncoding = ref_counts.get('operon_noncoding_nt', 0)
    
    add_line(f"   Protein sequences: {total_protein}")
    add_line(f"   Nucleotide gene sequences: {total_nt_genes}")
    add_line(f"   Non-coding regulatory sequences: {total_noncoding}")
    add_line(f"   Total query sequences: {total_protein + total_nt_genes + total_noncoding}")
    
    # BLAST search results
    add_line("\n2. BLAST Search Results by Strategy:")
    blast_stats = count_blast_files_and_hits()
    
    strategy_names = {
        "coding_protein": "Strategy A: tblastn (protein→nucleotide) vs Prokka genomes",
        "coding_nt": "Strategy B: blastn (nucleotide→nucleotide) vs Prokka genomes",
        "prokka_variants": "Strategy C: blastn vs Prokka-predicted CDS",
        "noncoding": "Strategy D: blastn for non-coding elements"
    }
    
    total_all_hits = 0
    for mode, stats in blast_stats.items():
        add_line(f"\n   {strategy_names.get(mode, mode)}:")
        add_line(f"     - Result files generated: {stats['files']}")
        add_line(f"     - Genomes with hits: {stats['genomes']:,}")
        add_line(f"     - Total BLAST hits: {stats['total_hits']:,}")
        total_all_hits += stats['total_hits']
    
    add_line(f"\n   Total hits across all strategies: {total_all_hits:,}")
    
    # Hit quality
    add_line("\n3. Hit Quality Analysis (Main Strategy):")
    quality_stats = analyze_hit_quality()
    if quality_stats and quality_stats.get('total_hits', 0) > 0:
        add_line(f"   Total hits analyzed: {quality_stats['total_hits']:,}")
        add_line(f"   High-quality hits (≥90% identity, ≥80% coverage): {quality_stats['high_quality_hits']:,}")
        add_line(f"   High-quality percentage: {quality_stats['high_quality_percent']:.1f}%")
        add_line(f"   Mean sequence identity: {quality_stats['mean_identity']:.1f}%")
        add_line(f"   Mean query coverage: {quality_stats['mean_coverage']:.1f}%")
    else:
        add_line("   No quality statistics available")
    
    # Operon completeness
    add_line("\n4. Operon Completeness Analysis:")
    completeness_stats = analyze_operon_completeness()
    if "total_genomes" in completeness_stats:
        add_line(f"   Total genomes analyzed: {completeness_stats['total_genomes']:,}")
        add_line(f"   Complete operons (7/7 genes): {completeness_stats['complete_operons']:,}")
        add_line(f"   Complete operon percentage: {completeness_stats['complete_percent']:.1f}%")
        if completeness_stats.get('has_promoter', 0) > 0:
            add_line(f"   Genomes with promoter detected: {completeness_stats['has_promoter']:,}")
        if completeness_stats.get('mean_completeness', 0) > 0:
            add_line(f"   Mean completeness score: {completeness_stats['mean_completeness']:.2f}")
    else:
        add_line("   Completeness analysis not available")
    
    # Target genomes count
    add_line("\n5. Target Genome Dataset:")
    add_line(f"   Total prokaryotic genomes searched: 8,287")
    add_line(f"   Genome annotation: Prokka-annotated E. faecalis genomes")
    add_line(f"   BLAST version: NCBI BLAST+ 2.12.0")
    
    # Summary for manuscript
    add_line("\n" + "=" * 60)
    add_line("MANUSCRIPT NUMBERS SUMMARY:")
    add_line("=" * 60)
    
    total_protein_queries = ref_counts.get('operon_genes_protein', 0)
    total_nt_queries = ref_counts.get('operon_genes_nt', 0) + ref_counts.get('operon_noncoding_nt', 0)
    total_genomes = blast_stats.get('coding_protein', {}).get('genomes', 0)
    
    add_line(f"Reference sequences: {total_protein_queries} protein, {total_nt_queries} nucleotide")
    add_line(f"Search strategies employed: 4 complementary approaches")
    add_line(f"Genomes searched: 8,287 prokaryotic genomes")
    add_line(f"Total BLAST result files: {sum(s.get('files', 0) for s in blast_stats.values()):,}")
    add_line(f"Total BLAST hits: {total_all_hits:,}")
    
    if quality_stats and quality_stats.get('total_hits', 0) > 0:
        add_line(f"High-quality hits (≥90% identity, ≥80% coverage): {quality_stats['high_quality_hits']:,} ({quality_stats['high_quality_percent']:.1f}%)")
        add_line(f"Mean sequence identity across all hits: {quality_stats['mean_identity']:.1f}%")
        add_line(f"Mean query coverage: {quality_stats['mean_coverage']:.1f}%")
    
    # Write to file if specified
    if output_file:
        try:
            with open(output_file, 'w') as f:
                for line in output_lines:
                    f.write(line + '\n')
            add_line(f"\nResults saved to: {output_file}")
        except Exception as e:
            add_line(f"\nError saving to file: {e}")
    
    return {
        'reference_counts': ref_counts,
        'blast_stats': blast_stats,
        'quality_stats': quality_stats,
        'completeness_stats': completeness_stats,
        'output_lines': output_lines
    }

if __name__ == "__main__":
    import sys
    
    # Change to script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    
    # Check for output file argument
    output_file = None
    if len(sys.argv) > 1:
        output_file = sys.argv[1]
    
    stats = generate_manuscript_stats(output_file)
