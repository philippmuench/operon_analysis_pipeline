#!/usr/bin/env python3
"""
Generate manuscript statistics from operon order analysis results.
"""

import os
import sys
import pandas as pd
import numpy as np
import json
import argparse

def generate_manuscript_stats(output_file=None):
    """Generate all statistics for manuscript.
    
    Args:
        output_file (str): Optional path to save output to file. If None, prints to console.
    """
    
    # Check if results exist
    results_file = "output/operon_order_analysis.csv"
    if not os.path.exists(results_file):
        print(f"Error: Results file not found: {results_file}")
        print("Please run analyze_operon_order.py first")
        return
    
    # Load results
    df = pd.read_csv(results_file)
    total_genomes = len(df)
    
    # Prepare output
    output_lines = []
    output_lines.append("="*60)
    output_lines.append("OPERON ORDER ANALYSIS - MANUSCRIPT STATISTICS")
    output_lines.append("="*60)
    output_lines.append("")
    
    # Dataset size
    output_lines.append("Dataset Overview:")
    output_lines.append(f"  Total genomes analyzed: {total_genomes:,}")
    output_lines.append("")
    
    # Operon completeness
    output_lines.append("Operon Completeness:")
    complete_operons = df['complete_operon'].sum()
    output_lines.append(f"  Complete operons (all 7 genes): {complete_operons:,} ({100*complete_operons/total_genomes:.1f}%)")
    
    # Gene count distribution
    gene_counts = df['n_genes_found'].value_counts().sort_index(ascending=False)
    output_lines.append(f"  Genomes by gene count:")
    for n_genes, count in gene_counts.items():
        percentage = 100 * count / total_genomes
        output_lines.append(f"    {n_genes} genes: {count:,} genomes ({percentage:.1f}%)")
    
    # Calculate genomes with near-complete operons
    near_complete = df[df['n_genes_found'] >= 5].shape[0]
    output_lines.append(f"  Near-complete operons (≥5 genes): {near_complete:,} ({100*near_complete/total_genomes:.1f}%)")
    output_lines.append("")
    
    # Individual gene presence
    output_lines.append("Individual Gene Presence:")
    gene_presence = {}
    canonical_order = ['frpC', 'glpC', 'ptsD', 'ptsC', 'ptsB', 'ptsA', 'fruR']
    
    for _, row in df.iterrows():
        genes = row['genes_found']
        if isinstance(genes, str):
            try:
                genes = eval(genes)
            except:
                genes = []
        for gene in genes:
            if gene in canonical_order:
                gene_presence[gene] = gene_presence.get(gene, 0) + 1
    
    for gene in canonical_order:
        count = gene_presence.get(gene, 0)
        percentage = 100 * count / total_genomes
        output_lines.append(f"  {gene}: {count:,} genomes ({percentage:.1f}%)")
    output_lines.append("")
    
    # Synteny conservation
    output_lines.append("Operon Organization:")
    
    # Same contig
    same_contig = df['same_contig'].sum()
    output_lines.append(f"  Genes on same contig: {same_contig:,} ({100*same_contig/total_genomes:.1f}%)")
    
    # Same strand
    same_strand = df['same_strand'].sum()
    output_lines.append(f"  Genes on same strand: {same_strand:,} ({100*same_strand/total_genomes:.1f}%)")
    
    # Canonical order
    canonical_order_count = df['canonical_order'].sum()
    output_lines.append(f"  Canonical gene order preserved: {canonical_order_count:,} ({100*canonical_order_count/total_genomes:.1f}%)")
    
    # Inversions and rearrangements
    inversions = df['has_inversion'].sum()
    rearrangements = df['has_rearrangement'].sum()
    output_lines.append(f"  With inversions: {inversions:,} ({100*inversions/total_genomes:.1f}%)")
    output_lines.append(f"  With rearrangements: {rearrangements:,} ({100*rearrangements/total_genomes:.1f}%)")
    output_lines.append("")
    
    # Synteny score statistics
    valid_synteny = df['synteny_score'].dropna()
    if len(valid_synteny) > 0:
        output_lines.append("Synteny Conservation Score:")
        output_lines.append(f"  Mean synteny score: {valid_synteny.mean():.3f}")
        output_lines.append(f"  Median synteny score: {valid_synteny.median():.3f}")
        output_lines.append(f"  High conservation (≥0.9): {(valid_synteny >= 0.9).sum():,} ({100*(valid_synteny >= 0.9).sum()/len(valid_synteny):.1f}%)")
        output_lines.append(f"  Moderate conservation (0.5-0.9): {((valid_synteny >= 0.5) & (valid_synteny < 0.9)).sum():,} ({100*((valid_synteny >= 0.5) & (valid_synteny < 0.9)).sum()/len(valid_synteny):.1f}%)")
        output_lines.append(f"  Low conservation (<0.5): {(valid_synteny < 0.5).sum():,} ({100*(valid_synteny < 0.5).sum()/len(valid_synteny):.1f}%)")
        output_lines.append("")
    
    # Operon type distribution
    output_lines.append("Operon Type Classification:")
    operon_types = df['operon_type'].value_counts()
    for operon_type, count in operon_types.items():
        percentage = 100 * count / total_genomes
        output_lines.append(f"  {operon_type}: {count:,} ({percentage:.1f}%)")
    output_lines.append("")
    
    # Inter-gene distances (for same-contig operons)
    valid_distances = df['mean_gene_distance'].dropna()
    if len(valid_distances) > 0:
        output_lines.append("Inter-gene Distances (same-contig operons):")
        output_lines.append(f"  Number of operons analyzed: {len(valid_distances):,}")
        output_lines.append(f"  Mean distance between genes: {valid_distances.mean():.0f} bp")
        output_lines.append(f"  Median distance: {valid_distances.median():.0f} bp")
        output_lines.append(f"  Maximum distance: {df['max_gene_distance'].max():.0f} bp")
        
        # Tight clustering
        tight_clusters = (valid_distances < 1000).sum()
        output_lines.append(f"  Tightly clustered (<1kb): {tight_clusters:,} ({100*tight_clusters/len(valid_distances):.1f}%)")
        output_lines.append("")
    
    # Number of contigs
    multi_contig = df[df['n_contigs'] > 1]
    if len(multi_contig) > 0:
        output_lines.append("Multi-contig Operons:")
        output_lines.append(f"  Operons spanning multiple contigs: {len(multi_contig):,} ({100*len(multi_contig)/total_genomes:.1f}%)")
        output_lines.append(f"  Mean number of contigs: {multi_contig['n_contigs'].mean():.1f}")
        output_lines.append(f"  Maximum contigs spanned: {multi_contig['n_contigs'].max()}")
        output_lines.append("")
    
    # Key findings for manuscript
    output_lines.append("Key Findings for Manuscript:")
    output_lines.append("-" * 40)
    
    # Most important statistics
    output_lines.append(f"• {complete_operons:,}/{total_genomes:,} genomes ({100*complete_operons/total_genomes:.1f}%) contain the complete 7-gene operon")
    
    if len(valid_synteny) > 0:
        high_synteny = (valid_synteny >= 0.9).sum()
        output_lines.append(f"• {high_synteny:,} genomes ({100*high_synteny/len(valid_synteny):.1f}%) maintain high synteny conservation (≥0.9)")
    
    output_lines.append(f"• {canonical_order_count:,} genomes ({100*canonical_order_count/total_genomes:.1f}%) preserve the canonical gene order")
    
    # PTS genes co-occurrence
    pts_genes = ['ptsA', 'ptsB', 'ptsC', 'ptsD']
    genomes_with_all_pts = 0
    for _, row in df.iterrows():
        genes = row['genes_found']
        if isinstance(genes, str):
            try:
                genes = eval(genes)
            except:
                genes = []
        if all(pts in genes for pts in pts_genes):
            genomes_with_all_pts += 1
    
    output_lines.append(f"• {genomes_with_all_pts:,} genomes ({100*genomes_with_all_pts/total_genomes:.1f}%) contain all four PTS components")
    
    # Complete and conserved
    complete_conserved = df[df['operon_type'] == 'Complete and conserved'].shape[0]
    output_lines.append(f"• {complete_conserved:,} genomes ({100*complete_conserved/total_genomes:.1f}%) have complete and fully conserved operons")
    
    output_lines.append("")
    output_lines.append("="*60)
    output_lines.append(f"Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Output results
    output_text = "\n".join(output_lines)
    
    if output_file:
        with open(output_file, 'w') as f:
            f.write(output_text)
        print(f"Manuscript statistics saved to: {output_file}")
    else:
        print(output_text)
    
    return output_text

def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="Generate manuscript statistics from operon order analysis")
    parser.add_argument("--output", "-o", default="output/manuscript_stats.txt",
                       help="Output file for statistics (default: output/manuscript_stats.txt)")
    
    args = parser.parse_args()
    
    # Generate statistics
    generate_manuscript_stats(args.output)

if __name__ == "__main__":
    main()