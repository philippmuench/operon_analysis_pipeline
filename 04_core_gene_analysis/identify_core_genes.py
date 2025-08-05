#!/usr/bin/env python3
"""
Identify core genes present in ≥95% of E. faecalis genomes.
"""

import os
import glob
from collections import Counter
from Bio import SeqIO
from multiprocessing import Pool, cpu_count
import pandas as pd
import argparse

def process_genome(prokka_dir):
    """Extract gene names from a single Prokka output."""
    genome_id = os.path.basename(prokka_dir)
    # Handle the .result naming convention
    gff_file = os.path.join(prokka_dir, f"{genome_id}.gff")
    
    genes = set()
    if os.path.exists(gff_file):
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 9 and parts[2] == 'gene':
                    # Extract gene name from attributes
                    attrs = parts[8]
                    for attr in attrs.split(';'):
                        if attr.startswith('gene='):
                            gene_name = attr.replace('gene=', '').strip()
                            genes.add(gene_name)
                            break
    
    return genome_id, genes

def main():
    # Setup paths
    output_dir = 'output'
    os.makedirs(output_dir, exist_ok=True)
    prokka_root = '../01_prokka_annotation/output/prokka_results'
    
    # Get all genome directories
    genome_dirs = glob.glob(os.path.join(prokka_root, '*/'))
    genome_dirs = [d.rstrip('/') for d in genome_dirs]  # Remove trailing slash
    genome_dirs.sort()  # Ensure consistent ordering
    
    n_genomes = len(genome_dirs)
    print(f"\n{'='*60}")
    print(f"Core Gene Identification")
    print(f"{'='*60}")
    print(f"Input directory: {prokka_root}")
    print(f"Number of genomes found: {n_genomes}")
    
    if n_genomes == 0:
        print("ERROR: No genomes found! Make sure to run Prokka annotation first.")
        return
    
    print(f"Processing mode: {'Test run (50 genomes)' if n_genomes <= 50 else f'Full run ({n_genomes} genomes)'}")
    print(f"{'='*60}\n")
    
    # Process genomes in parallel
    print("Extracting genes from all genomes...")
    with Pool(min(cpu_count() - 1, 40)) as pool:
        results = pool.map(process_genome, genome_dirs)
    
    # Count gene occurrences
    gene_counter = Counter()
    genome_genes = {}
    
    for genome_id, genes in results:
        genome_genes[genome_id] = genes
        gene_counter.update(genes)
    
    # Identify core genes (≥95% prevalence)
    threshold = 0.95
    min_genomes = int(n_genomes * threshold)
    
    core_genes = []
    gene_stats = []
    
    for gene, count in gene_counter.items():
        prevalence = count / n_genomes
        gene_stats.append({
            'gene': gene,
            'count': count,
            'prevalence': prevalence
        })
        
        if count >= min_genomes:
            core_genes.append(gene)
    
    # Save results
    # Core genes list
    core_genes_file = os.path.join(output_dir, 'core_genes_95pct.txt')
    with open(core_genes_file, 'w') as f:
        for gene in sorted(core_genes):
            f.write(f"{gene}\n")
    
    # Gene prevalence statistics
    stats_df = pd.DataFrame(gene_stats)
    if not stats_df.empty:
        stats_df = stats_df.sort_values('prevalence', ascending=False)
    stats_file = os.path.join(output_dir, 'gene_prevalence_stats.csv')
    stats_df.to_csv(stats_file, index=False)
    
    # Summary
    print(f"\nCore gene identification complete:")
    print(f"Total unique genes: {len(gene_counter)}")
    print(f"Core genes (≥{threshold*100}%): {len(core_genes)}")
    print(f"\nPrevalence distribution:")
    if not stats_df.empty:
        for threshold in [1.0, 0.99, 0.95, 0.90, 0.80, 0.50]:
            n = (stats_df['prevalence'] >= threshold).sum()
            print(f"  ≥{threshold*100:.0f}%: {n} genes")
    
    print(f"\nResults saved to:")
    print(f"  {core_genes_file}")
    print(f"  {stats_file}")

if __name__ == '__main__':
    main()