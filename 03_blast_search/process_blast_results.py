#!/usr/bin/env python3
"""
Process BLAST results to identify genomes with complete operons.
"""

import os
import pandas as pd
import glob
import argparse

def parse_blast_results(blast_dir, test_mode=False):
    """Parse all BLAST result files and compile hits."""
    all_results = []
    
    # Get all unique genome IDs from both types of BLAST files
    gene_files = glob.glob(os.path.join(blast_dir, '*_genes_blast.txt'))
    noncoding_files = glob.glob(os.path.join(blast_dir, '*_noncoding_blast.txt'))
    
    # Extract genome IDs
    genome_ids = set()
    for f in gene_files:
        genome_id = os.path.basename(f).replace('_genes_blast.txt', '')
        genome_ids.add(genome_id)
    for f in noncoding_files:
        genome_id = os.path.basename(f).replace('_noncoding_blast.txt', '')
        genome_ids.add(genome_id)
    
    genome_ids = sorted(list(genome_ids))
    
    if test_mode:
        genome_ids = genome_ids[:50]
        print(f"TEST MODE: Processing only first 50 genomes")
    
    print(f"Found {len(genome_ids)} genomes to process")
    
    for genome_id in genome_ids:
        # Process gene BLAST results
        gene_file = os.path.join(blast_dir, f'{genome_id}_genes_blast.txt')
        if os.path.exists(gene_file) and os.path.getsize(gene_file) > 0:
            df = pd.read_csv(gene_file, sep='\t', header=None,
                           names=['qseqid', 'sseqid', 'pident', 'length', 
                                 'mismatch', 'gapopen', 'qstart', 'qend',
                                 'sstart', 'send', 'evalue', 'bitscore', 'qcovs'])
            # Keep best hit per query gene
            best_hits = df.sort_values('bitscore', ascending=False).groupby('qseqid').first()
            best_hits['genome_id'] = genome_id
            best_hits['search_type'] = 'tblastn'
            all_results.append(best_hits)
        
        # Process non-coding BLAST results
        noncoding_file = os.path.join(blast_dir, f'{genome_id}_noncoding_blast.txt')
        if os.path.exists(noncoding_file) and os.path.getsize(noncoding_file) > 0:
            df = pd.read_csv(noncoding_file, sep='\t', header=None,
                           names=['qseqid', 'sseqid', 'pident', 'length', 
                                 'mismatch', 'gapopen', 'qstart', 'qend',
                                 'sstart', 'send', 'evalue', 'bitscore', 'qcovs'])
            # Keep best hit per query sequence
            best_hits = df.sort_values('bitscore', ascending=False).groupby('qseqid').first()
            best_hits['genome_id'] = genome_id
            best_hits['search_type'] = 'blastn'
            all_results.append(best_hits)
    
    return pd.concat(all_results) if all_results else pd.DataFrame()

def summarize_operon_presence(results_df, identity_threshold=90, coverage_threshold=80):
    """Summarize operon presence across genomes."""
    
    # Extract gene names from query IDs
    results_df['gene_name'] = results_df.index.str.split('|').str[1]
    
    # Define expected operon genes and non-coding elements
    operon_genes = [
        'frpC', 'glpC', 'ptsD', 'ptsC', 'ptsB', 'ptsA', 'fruR'
    ]
    noncoding_elements = ['promoter', 'pribnow_box']
    
    # Group by genome and check presence
    summary = []
    for genome_id, genome_data in results_df.groupby('genome_id'):
        genome_summary = {'genome_id': genome_id}
        
        # Check each operon gene
        genes_found = 0
        for gene in operon_genes:
            gene_hits = genome_data[genome_data['gene_name'] == gene]
            if not gene_hits.empty:
                best_hit = gene_hits.iloc[0]
                if best_hit['pident'] >= identity_threshold and best_hit['qcovs'] >= coverage_threshold:
                    genome_summary[f'{gene}_present'] = True
                    genome_summary[f'{gene}_identity'] = best_hit['pident']
                    genes_found += 1
                else:
                    genome_summary[f'{gene}_present'] = False
            else:
                genome_summary[f'{gene}_present'] = False
        
        genome_summary['operon_completeness'] = genes_found / len(operon_genes)
        genome_summary['n_genes_found'] = genes_found
        
        # Check non-coding elements
        for element in noncoding_elements:
            element_hits = genome_data[genome_data['gene_name'] == element]
            if not element_hits.empty:
                best_hit = element_hits.iloc[0]
                genome_summary[f'{element}_present'] = True
                genome_summary[f'{element}_identity'] = best_hit['pident']
            else:
                genome_summary[f'{element}_present'] = False
        
        summary.append(genome_summary)
    
    return pd.DataFrame(summary)

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Process BLAST results for operon analysis")
    args = parser.parse_args()
    
    # Paths: use in-repo structured locations
    blast_dir = 'output/blast_results'
    # Write summary to the core analysis output to keep results centralized
    output_dir = '../04_core_gene_analysis/output'
    os.makedirs(output_dir, exist_ok=True)
    
    # Parse BLAST results
    print("Parsing BLAST results...")
    results_df = parse_blast_results(blast_dir, test_mode=False)
    
    if results_df.empty:
        print("No BLAST results found!")
        return
    
    # Save raw results
    results_df.to_csv(os.path.join(output_dir, 'all_blast_hits.csv'))
    
    # Summarize operon presence
    print("Summarizing operon presence...")
    summary_df = summarize_operon_presence(results_df)
    
    # Save summary
    summary_df.to_csv(os.path.join(output_dir, 'operon_presence_summary.csv'), index=False)
    
    # Print statistics
    n_complete = (summary_df['operon_completeness'] == 1).sum()
    n_genomes = len(summary_df)
    
    print(f"\nOperon presence summary:")
    print(f"Total genomes analyzed: {n_genomes}")
    print(f"Genomes with complete operon: {n_complete} ({n_complete/n_genomes*100:.1f}%)")
    
    # Completeness distribution
    print("\nCompleteness distribution:")
    for completeness in [1.0, 0.85, 0.7, 0.5, 0]:
        n = (summary_df['operon_completeness'] >= completeness).sum()
        print(f"  ≥{completeness*100:.0f}% complete: {n} genomes ({n/n_genomes*100:.1f}%)")

if __name__ == '__main__':
    main()