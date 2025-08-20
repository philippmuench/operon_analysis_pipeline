#!/usr/bin/env python3
"""
Analyze operon gene order and synteny from BLAST mapping results.
Uses assembly-based BLAST results to determine gene positions and order.
"""

import os
import sys
import pandas as pd
import numpy as np
import glob
from collections import defaultdict, Counter
import json
import argparse
from pathlib import Path

# Define canonical operon gene order
CANONICAL_ORDER = ['frpC', 'glpC', 'ptsD', 'ptsC', 'ptsB', 'ptsA', 'fruR']
GENE_FUNCTIONS = {
    'frpC': 'Fructoselysine-6-phosphate deglycase',
    'glpC': 'Glucoselysine-6-phosphate deglycase',
    'ptsD': 'PTS system EIID component',
    'ptsC': 'PTS system EIIC component', 
    'ptsB': 'PTS system EIIB component',
    'ptsA': 'PTS system EIIA component',
    'fruR': 'Sigma-54 dependent transcriptional regulator'
}

def parse_blast_file(blast_file):
    """Parse a BLAST output file to extract hit information."""
    hits = []
    
    if not os.path.exists(blast_file):
        return hits
    
    try:
        # Read BLAST tabular output with 13 columns (includes qcovs)
        # Format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs
        df = pd.read_csv(blast_file, sep='\t', header=None,
                        names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                               'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                               'evalue', 'bitscore', 'qcovs'])
        
        for _, row in df.iterrows():
            hit = {
                'query': row['qseqid'],
                'subject': row['sseqid'],
                'identity': row['pident'],
                'start': min(row['sstart'], row['send']),
                'end': max(row['sstart'], row['send']),
                'strand': '+' if row['sstart'] < row['send'] else '-',
                'evalue': row['evalue'],
                'bitscore': row['bitscore']
            }
            hits.append(hit)
    
    except Exception as e:
        print(f"Error parsing {blast_file}: {e}")
    
    return hits

def extract_gene_from_query(query_id):
    """Extract gene name from query ID."""
    # Handle different query formats
    query_lower = query_id.lower()
    
    for gene in CANONICAL_ORDER:
        if gene.lower() in query_lower:
            return gene
    
    # Alternative names
    if 'fructoselysine' in query_lower and 'deglycase' in query_lower:
        return 'frpC'
    elif 'glucoselysine' in query_lower and 'deglycase' in query_lower:
        return 'glpC'
    elif 'pts' in query_lower:
        if 'eiid' in query_lower or 'iid' in query_lower:
            return 'ptsD'
        elif 'eiic' in query_lower or 'iic' in query_lower:
            return 'ptsC'
        elif 'eiib' in query_lower or 'iib' in query_lower:
            return 'ptsB'
        elif 'eiia' in query_lower or 'iia' in query_lower:
            return 'ptsA'
    elif 'sigma' in query_lower or 'regulator' in query_lower:
        return 'fruR'
    
    return None

def analyze_genome_operon_structure(genome_blast_dir):
    """Analyze operon structure for a single genome."""
    
    genome_id = os.path.basename(genome_blast_dir)
    
    # Find all BLAST result files for this genome
    blast_files = glob.glob(os.path.join(genome_blast_dir, '*.blast'))
    if not blast_files:
        blast_files = glob.glob(os.path.join(genome_blast_dir, '*.txt'))
    if not blast_files:
        blast_files = glob.glob(os.path.join(genome_blast_dir, '*'))
    
    # Collect all hits for operon genes
    gene_hits = defaultdict(list)
    
    for blast_file in blast_files:
        hits = parse_blast_file(blast_file)
        
        for hit in hits:
            gene = extract_gene_from_query(hit['query'])
            if gene and hit['identity'] >= 80 and hit['bitscore'] >= 50:  # Quality filters
                gene_hits[gene].append(hit)
    
    if not gene_hits:
        return {
            'genome_id': genome_id,
            'genes_found': [],
            'n_genes_found': 0,
            'complete_operon': False,
            'gene_order': [],
            'canonical_order': False,
            'same_contig': False,
            'same_strand': False,
            'n_contigs': 0,
            'synteny_score': 0,
            'has_inversion': False,
            'has_rearrangement': False,
            'mean_gene_distance': None,
            'max_gene_distance': None,
            'operon_type': 'No genes',
            'gene_details': None,
            'error': 'No operon genes found'
        }
    
    # Get best hit for each gene (highest bitscore)
    best_hits = {}
    for gene, hits in gene_hits.items():
        best_hit = max(hits, key=lambda x: x['bitscore'])
        best_hits[gene] = best_hit
    
    # Analyze operon structure
    genes_found = list(best_hits.keys())
    n_genes = len(genes_found)
    
    # Check completeness
    complete_operon = n_genes == len(CANONICAL_ORDER)
    
    # Analyze gene order if we have multiple genes
    gene_order = []
    gene_positions = []
    strands = []
    contigs = set()
    
    if n_genes > 0:
        # Sort genes by their genomic position
        gene_info = []
        for gene, hit in best_hits.items():
            gene_info.append({
                'gene': gene,
                'contig': hit['subject'],
                'start': hit['start'],
                'end': hit['end'],
                'strand': hit['strand'],
                'identity': hit['identity']
            })
            contigs.add(hit['subject'])
        
        # Sort by contig and position
        gene_info.sort(key=lambda x: (x['contig'], x['start']))
        
        gene_order = [g['gene'] for g in gene_info]
        gene_positions = [(g['start'], g['end']) for g in gene_info]
        strands = [g['strand'] for g in gene_info]
    
    # Check if genes are on same contig
    same_contig = len(contigs) == 1
    
    # Check if all genes are on same strand
    same_strand = len(set(strands)) == 1 if strands else False
    
    # Calculate synteny score (how well order matches canonical)
    synteny_score = calculate_synteny_score(gene_order, CANONICAL_ORDER)
    
    # Check for inversions or rearrangements
    has_inversion = False
    has_rearrangement = False
    
    if same_contig and n_genes > 1:
        # Check if order matches canonical (forward or reverse)
        canonical_indices = [CANONICAL_ORDER.index(g) for g in gene_order if g in CANONICAL_ORDER]
        if canonical_indices:
            # Check if indices are monotonic (either increasing or decreasing)
            diffs = np.diff(canonical_indices)
            if not all(d > 0 for d in diffs) and not all(d < 0 for d in diffs):
                has_rearrangement = True
        
        # Check for strand switches (inversions)
        if len(set(strands)) > 1:
            has_inversion = True
    
    # Calculate distances between adjacent genes
    gene_distances = []
    if same_contig and len(gene_positions) > 1:
        for i in range(len(gene_positions) - 1):
            dist = gene_positions[i+1][0] - gene_positions[i][1]
            gene_distances.append(dist)
    
    # Determine operon type
    if complete_operon and same_contig and same_strand and synteny_score > 0.9:
        operon_type = 'Complete and conserved'
    elif complete_operon and same_contig:
        operon_type = 'Complete with rearrangement'
    elif complete_operon:
        operon_type = 'Complete but fragmented'
    elif n_genes >= 5:
        operon_type = 'Near-complete'
    elif n_genes >= 3:
        operon_type = 'Partial'
    else:
        operon_type = 'Minimal'
    
    return {
        'genome_id': genome_id,
        'genes_found': genes_found,
        'n_genes_found': n_genes,
        'complete_operon': complete_operon,
        'gene_order': gene_order,
        'canonical_order': gene_order == CANONICAL_ORDER or gene_order == CANONICAL_ORDER[::-1],
        'same_contig': same_contig,
        'same_strand': same_strand,
        'n_contigs': len(contigs),
        'synteny_score': synteny_score,
        'has_inversion': has_inversion,
        'has_rearrangement': has_rearrangement,
        'mean_gene_distance': np.mean(gene_distances) if gene_distances else None,
        'max_gene_distance': max(gene_distances) if gene_distances else None,
        'operon_type': operon_type,
        'gene_details': json.dumps(gene_info) if 'gene_info' in locals() else None
    }

def analyze_genome_operon_structure_flat(blast_dir, genome_id):
    """Analyze operon structure for a single genome from flat BLAST directory."""
    
    # Find BLAST files for this genome
    # The genome_id already includes .result from the extraction
    genes_file = os.path.join(blast_dir, f"{genome_id}_genes_blast.txt")
    noncoding_file = os.path.join(blast_dir, f"{genome_id}_noncoding_blast.txt")
    
    # Collect all hits for operon genes
    gene_hits = defaultdict(list)
    
    # Parse genes BLAST file if it exists
    if os.path.exists(genes_file):
        hits = parse_blast_file(genes_file)
        
        for hit in hits:
            gene = extract_gene_from_query(hit['query'])
            if gene and hit['identity'] >= 80 and hit['bitscore'] >= 50:  # Quality filters
                gene_hits[gene].append(hit)
    
    if not gene_hits:
        return {
            'genome_id': genome_id,
            'genes_found': [],
            'n_genes_found': 0,
            'complete_operon': False,
            'gene_order': [],
            'canonical_order': False,
            'same_contig': False,
            'same_strand': False,
            'n_contigs': 0,
            'synteny_score': 0,
            'has_inversion': False,
            'has_rearrangement': False,
            'mean_gene_distance': None,
            'max_gene_distance': None,
            'operon_type': 'No genes',
            'gene_details': None,
            'error': 'No operon genes found'
        }
    
    # Get best hit for each gene (highest bitscore)
    best_hits = {}
    for gene, hits in gene_hits.items():
        best_hit = max(hits, key=lambda x: x['bitscore'])
        best_hits[gene] = best_hit
    
    # Analyze operon structure
    genes_found = list(best_hits.keys())
    n_genes = len(genes_found)
    
    # Check completeness
    complete_operon = n_genes == len(CANONICAL_ORDER)
    
    # Analyze gene order if we have multiple genes
    gene_order = []
    gene_positions = []
    strands = []
    contigs = set()
    
    if n_genes > 0:
        # Sort genes by their genomic position
        gene_info = []
        for gene, hit in best_hits.items():
            gene_info.append({
                'gene': gene,
                'contig': hit['subject'],
                'start': hit['start'],
                'end': hit['end'],
                'strand': hit['strand'],
                'identity': hit['identity']
            })
            contigs.add(hit['subject'])
        
        # Sort by contig and position
        gene_info.sort(key=lambda x: (x['contig'], x['start']))
        
        gene_order = [g['gene'] for g in gene_info]
        gene_positions = [(g['start'], g['end']) for g in gene_info]
        strands = [g['strand'] for g in gene_info]
    
    # Check if genes are on same contig
    same_contig = len(contigs) == 1
    
    # Check if all genes are on same strand
    same_strand = len(set(strands)) == 1 if strands else False
    
    # Calculate synteny score (how well order matches canonical)
    synteny_score = calculate_synteny_score(gene_order, CANONICAL_ORDER)
    
    # Check for inversions or rearrangements
    has_inversion = False
    has_rearrangement = False
    
    if same_contig and n_genes > 1:
        # Check if order matches canonical (forward or reverse)
        canonical_indices = [CANONICAL_ORDER.index(g) for g in gene_order if g in CANONICAL_ORDER]
        if canonical_indices:
            # Check if indices are monotonic (either increasing or decreasing)
            diffs = np.diff(canonical_indices)
            if not all(d > 0 for d in diffs) and not all(d < 0 for d in diffs):
                has_rearrangement = True
        
        # Check for strand switches (inversions)
        if len(set(strands)) > 1:
            has_inversion = True
    
    # Calculate distances between adjacent genes
    gene_distances = []
    if same_contig and len(gene_positions) > 1:
        for i in range(len(gene_positions) - 1):
            dist = gene_positions[i+1][0] - gene_positions[i][1]
            gene_distances.append(dist)
    
    # Determine operon type
    if complete_operon and same_contig and same_strand and synteny_score > 0.9:
        operon_type = 'Complete and conserved'
    elif complete_operon and same_contig:
        operon_type = 'Complete with rearrangement'
    elif complete_operon:
        operon_type = 'Complete but fragmented'
    elif n_genes >= 5:
        operon_type = 'Near-complete'
    elif n_genes >= 3:
        operon_type = 'Partial'
    else:
        operon_type = 'Minimal'
    
    return {
        'genome_id': genome_id,
        'genes_found': genes_found,
        'n_genes_found': n_genes,
        'complete_operon': complete_operon,
        'gene_order': gene_order,
        'canonical_order': gene_order == CANONICAL_ORDER or gene_order == CANONICAL_ORDER[::-1],
        'same_contig': same_contig,
        'same_strand': same_strand,
        'n_contigs': len(contigs),
        'synteny_score': synteny_score,
        'has_inversion': has_inversion,
        'has_rearrangement': has_rearrangement,
        'mean_gene_distance': np.mean(gene_distances) if gene_distances else None,
        'max_gene_distance': max(gene_distances) if gene_distances else None,
        'operon_type': operon_type,
        'gene_details': json.dumps(gene_info) if 'gene_info' in locals() else None
    }

def calculate_synteny_score(observed_order, canonical_order):
    """Calculate synteny conservation score."""
    if not observed_order:
        return 0
    
    # Get indices of observed genes in canonical order
    canonical_indices = []
    for gene in observed_order:
        if gene in canonical_order:
            canonical_indices.append(canonical_order.index(gene))
    
    if len(canonical_indices) < 2:
        return 1.0 if len(canonical_indices) == 1 else 0
    
    # Check if order is preserved (forward or reverse)
    forward_score = 0
    reverse_score = 0
    
    for i in range(len(canonical_indices) - 1):
        # Forward direction
        if canonical_indices[i+1] > canonical_indices[i]:
            forward_score += 1
        # Reverse direction  
        if canonical_indices[i+1] < canonical_indices[i]:
            reverse_score += 1
    
    max_score = max(forward_score, reverse_score)
    return max_score / (len(canonical_indices) - 1)

def main():
    parser = argparse.ArgumentParser(description="Analyze operon gene order from BLAST results")
    parser.add_argument("--blast-dir",
                       default="../03_blast_search/output/blast_results",
                       help="Directory containing BLAST results")
    parser.add_argument("--output-dir", default="output",
                       help="Output directory")
    parser.add_argument("--min-identity", type=float, default=80,
                       help="Minimum sequence identity threshold")
    parser.add_argument("--use-flat-structure", action="store_true", default=True,
                       help="BLAST results are in flat directory structure")
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("="*60)
    print("Operon Order and Synteny Analysis")
    print("="*60)
    print(f"BLAST results directory: {args.blast_dir}")
    print(f"Canonical operon order: {' -> '.join(CANONICAL_ORDER)}")
    print()
    
    # Handle different BLAST directory structures
    if args.use_flat_structure:
        # Files are directly in blast_dir with naming: ENT_*_genes_blast.txt
        blast_files = glob.glob(os.path.join(args.blast_dir, "ENT_*_genes_blast.txt"))
        
        if not blast_files:
            print(f"Error: No BLAST files found in {args.blast_dir}")
            print("Please run BLAST search (step 03) first")
            return 1
        
        # Extract unique genome IDs
        genome_ids = set()
        for bf in blast_files:
            basename = os.path.basename(bf)
            # Extract genome ID from filename (e.g., ENT_AA0002AA_AS.result from ENT_AA0002AA_AS.result_genes_blast.txt)
            genome_id = basename.replace('_genes_blast.txt', '')
            genome_ids.add(genome_id)
        
        print(f"Found {len(genome_ids)} genomes with BLAST results")
        
        # Process as flat structure
        results = []
        operon_types = Counter()
        
        for i, genome_id in enumerate(sorted(genome_ids)):
            if i % 500 == 0:
                print(f"  Processed {i}/{len(genome_ids)} genomes...")
            
            # Debug: Show first genome being processed
            if i == 0:
                print(f"  DEBUG: Processing genome ID: {genome_id}")
            
            # Create a pseudo-directory for compatibility
            result = analyze_genome_operon_structure_flat(args.blast_dir, genome_id)
            results.append(result)
            operon_types[result['operon_type']] += 1
            
            # Debug: Show first result
            if i == 0:
                print(f"  DEBUG: First result - genes found: {result.get('genes_found', [])})")
    
    else:
        # Original structure with subdirectories per genome
        genome_dirs = glob.glob(os.path.join(args.blast_dir, "ENT_*"))
        
        if not genome_dirs:
            print(f"Error: No genome directories found in {args.blast_dir}")
            print("Please run BLAST search (step 03) first")
            return 1
        
        print(f"Found {len(genome_dirs)} genome BLAST result directories")
        print("Analyzing operon structure...")
        
        # Analyze each genome
        results = []
        operon_types = Counter()
        
        for i, genome_dir in enumerate(genome_dirs):
            if i % 500 == 0:
                print(f"  Processed {i}/{len(genome_dirs)} genomes...")
            
            result = analyze_genome_operon_structure(genome_dir)
            results.append(result)
            operon_types[result['operon_type']] += 1
    
    # Create DataFrame
    df = pd.DataFrame(results)
    
    # Save detailed results
    output_file = os.path.join(args.output_dir, 'operon_order_analysis.csv')
    df.to_csv(output_file, index=False)
    print(f"\nDetailed results saved to: {output_file}")
    
    # Generate summary statistics
    print("\n" + "="*60)
    print("Summary Statistics")
    print("="*60)
    
    total_genomes = len(df)
    complete_operons = df['complete_operon'].sum()
    
    print(f"Total genomes analyzed: {total_genomes}")
    print(f"Genomes with complete operon (7 genes): {complete_operons} ({100*complete_operons/total_genomes:.1f}%)")
    print()
    
    # Gene presence statistics
    print("Gene presence across genomes:")
    gene_counts = Counter()
    for genes in df['genes_found']:
        for gene in genes:
            gene_counts[gene] += 1
    
    for gene in CANONICAL_ORDER:
        count = gene_counts.get(gene, 0)
        print(f"  {gene} ({GENE_FUNCTIONS[gene]}): {count} ({100*count/total_genomes:.1f}%)")
    
    print()
    
    # Operon organization statistics
    print("Operon organization:")
    same_contig = df['same_contig'].sum()
    same_strand = df['same_strand'].sum()
    canonical = df['canonical_order'].sum()
    
    print(f"  On same contig: {same_contig} ({100*same_contig/total_genomes:.1f}%)")
    print(f"  On same strand: {same_strand} ({100*same_strand/total_genomes:.1f}%)")
    print(f"  Canonical gene order: {canonical} ({100*canonical/total_genomes:.1f}%)")
    print(f"  With inversions: {df['has_inversion'].sum()} ({100*df['has_inversion'].sum()/total_genomes:.1f}%)")
    print(f"  With rearrangements: {df['has_rearrangement'].sum()} ({100*df['has_rearrangement'].sum()/total_genomes:.1f}%)")
    
    print()
    print("Operon types distribution:")
    for operon_type, count in sorted(operon_types.items(), key=lambda x: -x[1]):
        print(f"  {operon_type}: {count} ({100*count/total_genomes:.1f}%)")
    
    # Gene distance statistics
    valid_distances = df['mean_gene_distance'].dropna()
    if len(valid_distances) > 0:
        print()
        print("Inter-gene distances (for operons on same contig):")
        print(f"  Mean distance: {valid_distances.mean():.0f} bp")
        print(f"  Median distance: {valid_distances.median():.0f} bp")
        print(f"  Max distance: {df['max_gene_distance'].max():.0f} bp")
    
    # Save summary
    summary_file = os.path.join(args.output_dir, 'operon_order_summary.txt')
    with open(summary_file, 'w') as f:
        f.write("Operon Order and Synteny Analysis Summary\n")
        f.write("="*60 + "\n\n")
        f.write(f"Total genomes analyzed: {total_genomes}\n")
        f.write(f"Complete operons (7 genes): {complete_operons} ({100*complete_operons/total_genomes:.1f}%)\n\n")
        
        f.write("Gene presence:\n")
        for gene in CANONICAL_ORDER:
            count = gene_counts.get(gene, 0)
            f.write(f"  {gene}: {count} ({100*count/total_genomes:.1f}%)\n")
        
        f.write(f"\nOperon organization:\n")
        f.write(f"  Same contig: {same_contig} ({100*same_contig/total_genomes:.1f}%)\n")
        f.write(f"  Same strand: {same_strand} ({100*same_strand/total_genomes:.1f}%)\n")
        f.write(f"  Canonical order: {canonical} ({100*canonical/total_genomes:.1f}%)\n")
        
        f.write(f"\nOperon types:\n")
        for operon_type, count in sorted(operon_types.items(), key=lambda x: -x[1]):
            f.write(f"  {operon_type}: {count} ({100*count/total_genomes:.1f}%)\n")
    
    print(f"\nSummary saved to: {summary_file}")
    
    # Create gene co-occurrence matrix
    print("\nCalculating gene co-occurrence...")
    cooccurrence = pd.DataFrame(0, index=CANONICAL_ORDER, columns=CANONICAL_ORDER)
    
    for genes in df['genes_found']:
        for i, gene1 in enumerate(genes):
            for gene2 in genes:
                if gene1 in CANONICAL_ORDER and gene2 in CANONICAL_ORDER:
                    cooccurrence.loc[gene1, gene2] += 1
    
    cooccurrence_file = os.path.join(args.output_dir, 'gene_cooccurrence_matrix.csv')
    cooccurrence.to_csv(cooccurrence_file)
    print(f"Gene co-occurrence matrix saved to: {cooccurrence_file}")
    
    print("\n" + "="*60)
    print("Analysis complete!")

if __name__ == "__main__":
    main()