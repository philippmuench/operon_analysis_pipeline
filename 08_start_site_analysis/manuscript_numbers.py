#!/usr/bin/env python3
"""
Generate manuscript statistics from start site analysis results.
"""

import os
import sys
import pandas as pd
import numpy as np
import argparse
from collections import Counter

def generate_manuscript_stats(output_file=None):
    """Generate all statistics for manuscript.
    
    Args:
        output_file (str): Optional path to save output to file. If None, prints to console.
    """
    
    # Check if results exist
    summary_file = "output/start_site_summary.tsv"
    if not os.path.exists(summary_file):
        print(f"Error: Results file not found: {summary_file}")
        print("Please run analyze_start_sites.py first")
        return
    
    # Load results
    df = pd.read_csv(summary_file, sep='\t')
    
    # Prepare output
    output_lines = []
    output_lines.append("="*60)
    output_lines.append("START SITE ANALYSIS - MANUSCRIPT STATISTICS")
    output_lines.append("="*60)
    output_lines.append("")
    
    # Dataset overview
    output_lines.append("Dataset Overview:")
    total_genes = len(df)
    unique_genomes = df['genome_id'].nunique()
    unique_gene_types = df['gene'].nunique()
    output_lines.append(f"  Total gene instances analyzed: {total_genes:,}")
    output_lines.append(f"  Unique genomes: {unique_genomes:,}")
    output_lines.append(f"  Unique operon genes analyzed: {unique_gene_types}")
    output_lines.append("")
    
    # Start codon usage
    output_lines.append("Start Codon Usage:")
    start_codon_counts = df['start_codon'].value_counts()
    total_valid_starts = start_codon_counts.sum()
    
    for codon, count in start_codon_counts.items():
        percentage = 100 * count / total_valid_starts
        output_lines.append(f"  {codon}: {count:,} ({percentage:.1f}%)")
    output_lines.append("")
    
    # Alternative start sites
    output_lines.append("Alternative Start Sites:")
    upstream_found = df['upstream_candidate_found'].sum()
    upstream_pct = 100 * upstream_found / total_genes
    output_lines.append(f"  Genes with upstream alternative starts: {upstream_found:,} ({upstream_pct:.1f}%)")
    
    # RBS analysis for upstream starts
    upstream_df = df[df['upstream_candidate_found'] == True]
    if len(upstream_df) > 0:
        upstream_with_rbs = upstream_df['upstream_has_rbs'].sum()
        upstream_rbs_pct = 100 * upstream_with_rbs / len(upstream_df)
        output_lines.append(f"  Upstream starts with RBS: {upstream_with_rbs:,} ({upstream_rbs_pct:.1f}%)")
    
    # Average offset when upstream start exists
    if 'upstream_offset_nt' in df.columns:
        valid_offsets = df['upstream_offset_nt'].dropna()
        if len(valid_offsets) > 0:
            output_lines.append(f"  Mean upstream offset: {valid_offsets.mean():.1f} nt")
            output_lines.append(f"  Median upstream offset: {valid_offsets.median():.0f} nt")
            output_lines.append(f"  Max upstream offset: {valid_offsets.max():.0f} nt")
    output_lines.append("")
    
    # Classification breakdown
    output_lines.append("Start Site Classifications:")
    if 'classification' in df.columns:
        class_counts = df['classification'].value_counts()
        for classification, count in class_counts.items():
            percentage = 100 * count / total_genes
            output_lines.append(f"  {classification}: {count:,} ({percentage:.1f}%)")
    output_lines.append("")
    
    # RBS motif analysis
    output_lines.append("Ribosome Binding Site (RBS) Analysis:")
    
    # RBS presence from GFF annotations
    if 'gff_rbs_motif' in df.columns:
        has_gff_rbs = df['gff_rbs_motif'].notna().sum()
        gff_rbs_pct = 100 * has_gff_rbs / total_genes
        output_lines.append(f"  Genes with annotated RBS: {has_gff_rbs:,} ({gff_rbs_pct:.1f}%)")
        
        # RBS motif types
        rbs_motifs = df['gff_rbs_motif'].dropna().value_counts()
        if len(rbs_motifs) > 0:
            output_lines.append("  Common RBS motifs:")
            for motif, count in rbs_motifs.head(5).items():
                output_lines.append(f"    {motif}: {count:,}")
    
    # RBS spacer distances
    if 'gff_rbs_spacer' in df.columns:
        valid_spacers = pd.to_numeric(df['gff_rbs_spacer'], errors='coerce').dropna()
        if len(valid_spacers) > 0:
            output_lines.append(f"  Mean RBS-start spacer: {valid_spacers.mean():.1f} nt")
            output_lines.append(f"  Typical range: {valid_spacers.quantile(0.25):.0f}-{valid_spacers.quantile(0.75):.0f} nt")
    output_lines.append("")
    
    # Gene-specific analysis
    output_lines.append("Gene-Specific Patterns:")
    gene_stats = []
    for gene in df['gene'].unique():
        gene_df = df[df['gene'] == gene]
        gene_total = len(gene_df)
        
        # Start codon preference
        most_common_start = gene_df['start_codon'].value_counts().index[0] if len(gene_df) > 0 else 'NA'
        atg_count = (gene_df['start_codon'] == 'ATG').sum()
        atg_pct = 100 * atg_count / gene_total if gene_total > 0 else 0
        
        # Alternative starts
        alt_starts = gene_df['upstream_candidate_found'].sum()
        alt_pct = 100 * alt_starts / gene_total if gene_total > 0 else 0
        
        gene_stats.append({
            'gene': gene,
            'total': gene_total,
            'primary_start': most_common_start,
            'atg_pct': atg_pct,
            'alt_start_pct': alt_pct
        })
    
    gene_stats_df = pd.DataFrame(gene_stats).sort_values('gene')
    for _, row in gene_stats_df.iterrows():
        output_lines.append(f"  {row['gene']}: {row['primary_start']} start ({row['atg_pct']:.0f}% ATG), {row['alt_start_pct']:.0f}% have alt starts")
    output_lines.append("")
    
    # Partial genes
    if 'gff_partial' in df.columns:
        partial_genes = df['gff_partial'].notna().sum()
        if partial_genes > 0:
            partial_pct = 100 * partial_genes / total_genes
            output_lines.append("Partial/Truncated Genes:")
            output_lines.append(f"  Partial gene predictions: {partial_genes:,} ({partial_pct:.1f}%)")
            output_lines.append("")
    
    # Key findings for manuscript
    output_lines.append("Key Findings for Manuscript:")
    output_lines.append("-" * 40)
    
    # Most important statistics
    atg_total = (df['start_codon'] == 'ATG').sum()
    atg_percentage = 100 * atg_total / total_valid_starts if total_valid_starts > 0 else 0
    output_lines.append(f"• {atg_total:,}/{total_valid_starts:,} ({atg_percentage:.1f}%) genes use ATG as start codon")
    
    gtg_total = (df['start_codon'] == 'GTG').sum()
    gtg_percentage = 100 * gtg_total / total_valid_starts if total_valid_starts > 0 else 0
    ttg_total = (df['start_codon'] == 'TTG').sum()
    ttg_percentage = 100 * ttg_total / total_valid_starts if total_valid_starts > 0 else 0
    
    if gtg_total > 0 or ttg_total > 0:
        output_lines.append(f"• Alternative starts: GTG ({gtg_percentage:.1f}%), TTG ({ttg_percentage:.1f}%)")
    
    output_lines.append(f"• {upstream_found:,} ({upstream_pct:.1f}%) genes have potential upstream start sites")
    
    if 'upstream_has_rbs' in df.columns and upstream_found > 0:
        upstream_with_rbs = df[df['upstream_candidate_found'] == True]['upstream_has_rbs'].sum()
        output_lines.append(f"• {upstream_with_rbs:,} upstream sites have associated RBS motifs")
    
    # Classification insights
    if 'classification' in df.columns:
        alt_upstream_rbs = (df['classification'] == 'Alt_upstream_with_RBS').sum()
        if alt_upstream_rbs > 0:
            alt_upstream_rbs_pct = 100 * alt_upstream_rbs / total_genes
            output_lines.append(f"• {alt_upstream_rbs:,} ({alt_upstream_rbs_pct:.1f}%) genes likely use alternative upstream starts with RBS")
    
    # PTS genes specifically (often have alternative starts)
    pts_genes = ['ptsA', 'ptsB', 'ptsC', 'ptsD']
    pts_df = df[df['gene'].isin(pts_genes)]
    if len(pts_df) > 0:
        pts_alt_starts = pts_df['upstream_candidate_found'].sum()
        pts_alt_pct = 100 * pts_alt_starts / len(pts_df)
        output_lines.append(f"• PTS genes show {pts_alt_pct:.1f}% frequency of alternative start sites")
    
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
    parser = argparse.ArgumentParser(description="Generate manuscript statistics from start site analysis")
    parser.add_argument("--output", "-o", default="output/manuscript_stats.txt",
                       help="Output file for statistics (default: output/manuscript_stats.txt)")
    
    args = parser.parse_args()
    
    # Generate statistics
    generate_manuscript_stats(args.output)

if __name__ == "__main__":
    main()