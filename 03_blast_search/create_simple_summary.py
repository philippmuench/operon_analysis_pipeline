#!/usr/bin/env python3
"""
Create a simple summary table showing operon presence/absence per genome.
"""

import os
import pandas as pd
import glob
import argparse

def create_simple_summary(blast_dir, output_dir, test_mode=False):
    """Create a simple presence/absence matrix for operon components."""
    
    # Get all unique genome IDs
    gene_files = glob.glob(os.path.join(blast_dir, '*_genes_blast.txt'))
    genome_ids = []
    for f in gene_files:
        genome_id = os.path.basename(f).replace('_genes_blast.txt', '')
        genome_ids.append(genome_id)
    
    genome_ids.sort()
    
    if test_mode:
        genome_ids = genome_ids[:50]
    
    print(f"Processing {len(genome_ids)} genomes...")
    
    # Define components
    gene_components = ['frpC', 'glpC', 'ptsD', 'ptsC', 'ptsB', 'ptsA', 'fruR']
    noncoding_components = ['promoter']  # Skip pribnow_box as it's too short
    all_components = gene_components + noncoding_components
    
    # Initialize results
    summary_data = []
    
    for genome_id in genome_ids:
        genome_data = {'genome_id': genome_id}
        
        # Process gene BLAST results
        gene_file = os.path.join(blast_dir, f'{genome_id}_genes_blast.txt')
        if os.path.exists(gene_file) and os.path.getsize(gene_file) > 0:
            df = pd.read_csv(gene_file, sep='\t', header=None,
                           names=['qseqid', 'sseqid', 'pident', 'length', 
                                 'mismatch', 'gapopen', 'qstart', 'qend',
                                 'sstart', 'send', 'evalue', 'bitscore', 'qcovs'])
            
            # Extract gene names and get best hit per gene
            df['gene_name'] = df['qseqid'].str.split('|').str[1]
            best_hits = df.sort_values('bitscore', ascending=False).groupby('gene_name').first()
            
            for gene in gene_components:
                if gene in best_hits.index:
                    hit = best_hits.loc[gene]
                    if hit['pident'] >= 90 and hit['qcovs'] >= 80:
                        genome_data[gene] = 1
                        genome_data[f'{gene}_identity'] = hit['pident']
                        genome_data[f'{gene}_coverage'] = hit['qcovs']
                    else:
                        genome_data[gene] = 0
                else:
                    genome_data[gene] = 0
        
        # Process non-coding BLAST results
        noncoding_file = os.path.join(blast_dir, f'{genome_id}_noncoding_blast.txt')
        if os.path.exists(noncoding_file) and os.path.getsize(noncoding_file) > 0:
            df = pd.read_csv(noncoding_file, sep='\t', header=None,
                           names=['qseqid', 'sseqid', 'pident', 'length', 
                                 'mismatch', 'gapopen', 'qstart', 'qend',
                                 'sstart', 'send', 'evalue', 'bitscore', 'qcovs'])
            
            df['element_name'] = df['qseqid'].str.split('|').str[1]
            best_hits = df.sort_values('bitscore', ascending=False).groupby('element_name').first()
            
            for element in noncoding_components:
                if element in best_hits.index:
                    hit = best_hits.loc[element]
                    genome_data[element] = 1
                    genome_data[f'{element}_identity'] = hit['pident']
                    genome_data[f'{element}_coverage'] = hit['qcovs']
                else:
                    genome_data[element] = 0
        
        # Calculate completeness
        genes_present = sum(genome_data.get(gene, 0) for gene in gene_components)
        genome_data['n_genes'] = genes_present
        genome_data['gene_completeness'] = genes_present / len(gene_components) * 100
        genome_data['has_promoter'] = genome_data.get('promoter', 0)
        genome_data['is_complete'] = 1 if genes_present == len(gene_components) else 0
        
        summary_data.append(genome_data)
    
    # Create DataFrame
    summary_df = pd.DataFrame(summary_data)
    
    # Reorder columns
    basic_cols = ['genome_id', 'n_genes', 'gene_completeness', 'is_complete', 'has_promoter']
    component_cols = all_components
    detail_cols = []
    for comp in all_components:
        if f'{comp}_identity' in summary_df.columns:
            detail_cols.extend([f'{comp}_identity', f'{comp}_coverage'])
    
    ordered_cols = basic_cols + component_cols + detail_cols
    summary_df = summary_df[ordered_cols]
    
    # Save summary
    output_file = os.path.join(output_dir, 'operon_simple_summary.csv')
    summary_df.to_csv(output_file, index=False)
    
    # Print statistics
    print(f"\n=== OPERON PRESENCE SUMMARY ===")
    print(f"Total genomes: {len(summary_df)}")
    print(f"Complete operons (all 7 genes): {summary_df['is_complete'].sum()} ({summary_df['is_complete'].sum()/len(summary_df)*100:.1f}%)")
    print(f"Genomes with promoter: {summary_df['has_promoter'].sum()} ({summary_df['has_promoter'].sum()/len(summary_df)*100:.1f}%)")
    
    print(f"\nGene prevalence:")
    for gene in gene_components:
        count = summary_df[gene].sum()
        print(f"  {gene:8s}: {count:3d} ({count/len(summary_df)*100:5.1f}%)")
    
    print(f"\nCompleteness distribution:")
    for n in range(8):
        count = (summary_df['n_genes'] == n).sum()
        if count > 0:
            print(f"  {n} genes: {count:3d} genomes ({count/len(summary_df)*100:5.1f}%)")
    
    print(f"\nSummary saved to: {output_file}")
    
    # Also create a presence/absence matrix
    pa_matrix = summary_df[['genome_id'] + all_components].copy()
    pa_file = os.path.join(output_dir, 'operon_presence_absence_matrix.csv')
    pa_matrix.to_csv(pa_file, index=False)
    print(f"Presence/absence matrix saved to: {pa_file}")
    
    return summary_df

def main():
    parser = argparse.ArgumentParser(description="Create simple operon summary")
    args = parser.parse_args()
    
    # Paths: read local BLAST results, write to core analysis outputs
    blast_dir = 'output/blast_results'
    output_dir = '../04_core_gene_analysis/output'
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Create summary
    summary_df = create_simple_summary(blast_dir, output_dir, test_mode=False)

if __name__ == '__main__':
    main()