#!/usr/bin/env python3
"""
Create overview of BLAST results showing prevalence of each operon component.
"""

import os
import pandas as pd
import glob
import argparse

def parse_all_blast_hits(blast_dir, test_mode=False):
    """Parse all BLAST result files to get all hits (not just best)."""
    all_hits = []
    
    # Get all unique genome IDs
    gene_files = glob.glob(os.path.join(blast_dir, '*_genes_blast.txt'))
    noncoding_files = glob.glob(os.path.join(blast_dir, '*_noncoding_blast.txt'))
    
    genome_ids = set()
    for f in gene_files:
        genome_id = os.path.basename(f).replace('_genes_blast.txt', '')
        genome_ids.add(genome_id)
    
    genome_ids = sorted(list(genome_ids))
    
    if test_mode:
        genome_ids = genome_ids[:50]
        print(f"TEST MODE: Processing only first 50 genomes")
    
    print(f"Processing {len(genome_ids)} genomes...")
    
    for genome_id in genome_ids:
        # Process gene BLAST results
        gene_file = os.path.join(blast_dir, f'{genome_id}_genes_blast.txt')
        if os.path.exists(gene_file) and os.path.getsize(gene_file) > 0:
            df = pd.read_csv(gene_file, sep='\t', header=None,
                           names=['qseqid', 'sseqid', 'pident', 'length', 
                                 'mismatch', 'gapopen', 'qstart', 'qend',
                                 'sstart', 'send', 'evalue', 'bitscore', 'qcovs'])
            df['genome_id'] = genome_id
            df['search_type'] = 'tblastn'
            all_hits.append(df)
        
        # Process non-coding BLAST results
        noncoding_file = os.path.join(blast_dir, f'{genome_id}_noncoding_blast.txt')
        if os.path.exists(noncoding_file) and os.path.getsize(noncoding_file) > 0:
            df = pd.read_csv(noncoding_file, sep='\t', header=None,
                           names=['qseqid', 'sseqid', 'pident', 'length', 
                                 'mismatch', 'gapopen', 'qstart', 'qend',
                                 'sstart', 'send', 'evalue', 'bitscore', 'qcovs'])
            df['genome_id'] = genome_id
            df['search_type'] = 'blastn'
            all_hits.append(df)
    
    return pd.concat(all_hits) if all_hits else pd.DataFrame()

def create_overview(results_df, output_dir):
    """Create overview statistics and visualizations."""
    
    # Extract gene/element names from query IDs
    results_df['element_name'] = results_df['qseqid'].str.split('|').str[1]
    
    # Define all operon components
    operon_components = {
        'Genes': ['frpC', 'glpC', 'ptsD', 'ptsC', 'ptsB', 'ptsA', 'fruR'],
        'Non-coding': ['promoter', 'pribnow_box']
    }
    all_components = operon_components['Genes'] + operon_components['Non-coding']
    
    # 1. Overall statistics
    print("\n=== OVERALL STATISTICS ===")
    print(f"Total genomes analyzed: {results_df['genome_id'].nunique()}")
    print(f"Total BLAST hits: {len(results_df)}")
    
    # 2. Hits per genome
    hits_per_genome = results_df.groupby('genome_id').size()
    print(f"\nHits per genome:")
    print(f"  Mean: {hits_per_genome.mean():.1f}")
    print(f"  Median: {hits_per_genome.median():.0f}")
    print(f"  Range: {hits_per_genome.min()}-{hits_per_genome.max()}")
    
    # 3. Component prevalence (using high quality hits)
    print("\n=== COMPONENT PREVALENCE ===")
    print("(Using hits with ≥90% identity and ≥80% coverage)\n")
    
    prevalence_data = []
    for component in all_components:
        component_hits = results_df[results_df['element_name'] == component]
        high_quality = component_hits[(component_hits['pident'] >= 90) & 
                                     (component_hits['qcovs'] >= 80)]
        genomes_with_component = high_quality['genome_id'].nunique()
        total_genomes = results_df['genome_id'].nunique()
        prevalence = genomes_with_component / total_genomes * 100
        
        prevalence_data.append({
            'Component': component,
            'Type': 'Gene' if component in operon_components['Genes'] else 'Non-coding',
            'Genomes_with_component': genomes_with_component,
            'Total_genomes': total_genomes,
            'Prevalence_%': prevalence,
            'Total_hits': len(component_hits),
            'High_quality_hits': len(high_quality)
        })
        
        print(f"{component:15s}: {prevalence:6.2f}% ({genomes_with_component}/{total_genomes} genomes)")
    
    prevalence_df = pd.DataFrame(prevalence_data)
    prevalence_df.to_csv(os.path.join(output_dir, 'component_prevalence.csv'), index=False)
    
    # 4. Multiple hits analysis
    print("\n=== MULTIPLE HITS ANALYSIS ===")
    multiple_hits = results_df.groupby(['genome_id', 'element_name']).size()
    for component in all_components:
        if component in multiple_hits.index.get_level_values('element_name'):
            component_counts = multiple_hits[:, component]
            genomes_with_multiple = (component_counts > 1).sum()
            if genomes_with_multiple > 0:
                max_copies = component_counts.max()
                print(f"{component}: {genomes_with_multiple} genomes have multiple hits (max: {max_copies})")
    
    # 5. Identity distribution per component
    print("\n=== IDENTITY DISTRIBUTION ===")
    identity_stats = []
    for component in all_components:
        component_hits = results_df[results_df['element_name'] == component]
        if not component_hits.empty:
            stats = {
                'Component': component,
                'Mean_identity': component_hits['pident'].mean(),
                'Median_identity': component_hits['pident'].median(),
                'Min_identity': component_hits['pident'].min(),
                'Max_identity': component_hits['pident'].max()
            }
            identity_stats.append(stats)
            print(f"{component:15s}: mean={stats['Mean_identity']:.1f}%, median={stats['Median_identity']:.1f}%")
    
    identity_df = pd.DataFrame(identity_stats)
    identity_df.to_csv(os.path.join(output_dir, 'identity_statistics.csv'), index=False)
    
    return prevalence_df

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Create overview of BLAST results")
    args = parser.parse_args()
    
    # Paths: use local BLAST results, write to core analysis outputs
    blast_dir = 'output/blast_results'
    output_dir = '../04_core_gene_analysis/output'
    
    # Parse all BLAST hits
    print("Parsing BLAST results...")
    results_df = parse_all_blast_hits(blast_dir, test_mode=False)
    
    if results_df.empty:
        print("No BLAST results found!")
        return
    
    # Create overview
    prevalence_df = create_overview(results_df, output_dir)
    
    # Save complete results
    results_df.to_csv(os.path.join(output_dir, 'all_blast_hits_complete.csv'), index=False)
    print(f"\nComplete results saved to: {os.path.join(output_dir, 'all_blast_hits_complete.csv')}")

if __name__ == '__main__':
    main()