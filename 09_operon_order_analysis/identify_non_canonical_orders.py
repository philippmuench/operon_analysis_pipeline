#!/usr/bin/env python3
"""
Identify genomes with non-canonical operon gene order.
Lists all cases where the gene order differs from the reference.
"""

import pandas as pd
import json
import sys

# Define canonical order
CANONICAL_ORDER = ['frpC', 'glpC', 'ptsD', 'ptsC', 'ptsB', 'ptsA', 'fruR']

def main():
    # Load the analysis results
    results_file = "output/operon_order_analysis.csv"
    
    print("Loading operon order analysis results...")
    df = pd.read_csv(results_file)
    
    # Filter for genomes with complete operons
    complete_operons = df[df['complete_operon'] == True]
    print(f"Found {len(complete_operons)} genomes with complete operons")
    
    # Find non-canonical orders
    non_canonical = []
    canonical_count = 0
    
    for idx, row in complete_operons.iterrows():
        genome_id = row['genome_id']
        gene_order = row['gene_order']
        
        # Parse gene order (it might be stored as string representation of list)
        if isinstance(gene_order, str):
            try:
                gene_order = eval(gene_order)
            except:
                print(f"Warning: Could not parse gene order for {genome_id}: {gene_order}")
                continue
        
        # Check if it matches canonical order (forward or reverse)
        if gene_order == CANONICAL_ORDER:
            canonical_count += 1
        elif gene_order == CANONICAL_ORDER[::-1]:
            canonical_count += 1
            # Note: reverse order is still considered canonical (opposite strand)
        else:
            non_canonical.append({
                'genome_id': genome_id,
                'gene_order': gene_order,
                'n_genes': len(gene_order),
                'same_contig': row.get('same_contig', None),
                'same_strand': row.get('same_strand', None),
                'synteny_score': row.get('synteny_score', None)
            })
    
    print(f"\nCanonical order (forward or reverse): {canonical_count}")
    print(f"Non-canonical order: {len(non_canonical)}")
    
    # Also check for incomplete operons with unusual patterns
    partial_operons = df[(df['n_genes_found'] > 0) & (df['n_genes_found'] < 7)]
    print(f"\nPartial operons (1-6 genes): {len(partial_operons)}")
    
    # Write results to file
    output_file = "output/non_canonical_gene_orders.txt"
    with open(output_file, 'w') as f:
        f.write("="*70 + "\n")
        f.write("NON-CANONICAL OPERON GENE ORDERS\n")
        f.write("="*70 + "\n\n")
        
        f.write(f"Canonical order: {' -> '.join(CANONICAL_ORDER)}\n")
        f.write(f"Reverse canonical: {' -> '.join(CANONICAL_ORDER[::-1])}\n\n")
        
        f.write(f"Total genomes analyzed: {len(df)}\n")
        f.write(f"Complete operons: {len(complete_operons)}\n")
        f.write(f"Canonical order (forward or reverse): {canonical_count}\n")
        f.write(f"Non-canonical order: {len(non_canonical)}\n\n")
        
        if non_canonical:
            f.write("-"*70 + "\n")
            f.write("GENOMES WITH NON-CANONICAL GENE ORDER:\n")
            f.write("-"*70 + "\n\n")
            
            for i, case in enumerate(non_canonical, 1):
                f.write(f"Case {i}:\n")
                f.write(f"  Genome: {case['genome_id']}\n")
                f.write(f"  Gene order: {' -> '.join(case['gene_order'])}\n")
                f.write(f"  Number of genes: {case['n_genes']}\n")
                f.write(f"  Same contig: {case['same_contig']}\n")
                f.write(f"  Same strand: {case['same_strand']}\n")
                f.write(f"  Synteny score: {case['synteny_score']:.3f}\n")
                
                # Identify what's different
                differences = []
                for j, gene in enumerate(case['gene_order']):
                    if j < len(CANONICAL_ORDER):
                        if gene != CANONICAL_ORDER[j]:
                            differences.append(f"Position {j+1}: expected {CANONICAL_ORDER[j]}, found {gene}")
                
                if differences:
                    f.write(f"  Differences from canonical:\n")
                    for diff in differences:
                        f.write(f"    - {diff}\n")
                
                f.write("\n")
        else:
            f.write("All complete operons have canonical gene order!\n")
        
        # Add section for partial operons
        if len(partial_operons) > 0:
            f.write("\n" + "-"*70 + "\n")
            f.write("PARTIAL OPERONS:\n")
            f.write("-"*70 + "\n\n")
            
            for idx, row in partial_operons.iterrows():
                genome_id = row['genome_id']
                genes = row['genes_found']
                if isinstance(genes, str):
                    try:
                        genes = eval(genes)
                    except:
                        genes = []
                
                gene_order = row['gene_order']
                if isinstance(gene_order, str):
                    try:
                        gene_order = eval(gene_order)
                    except:
                        gene_order = []
                
                f.write(f"Genome: {genome_id}\n")
                f.write(f"  Genes present ({row['n_genes_found']}): {', '.join(genes)}\n")
                if gene_order:
                    f.write(f"  Gene order: {' -> '.join(gene_order)}\n")
                f.write("\n")
    
    print(f"\nResults written to: {output_file}")
    
    # Additional check for potential errors in our analysis
    print("\n" + "="*70)
    print("QUALITY CHECK:")
    print("="*70)
    
    # Check if canonical_order field matches our manual check
    canonical_field = complete_operons[complete_operons['canonical_order'] == True]
    print(f"Genomes marked as canonical in data: {len(canonical_field)}")
    print(f"Genomes we counted as canonical: {canonical_count}")
    
    if len(canonical_field) != canonical_count:
        print("WARNING: Discrepancy detected! Investigating...")
        
        # Find discrepancies
        for idx, row in complete_operons.iterrows():
            genome_id = row['genome_id']
            gene_order = row['gene_order']
            marked_canonical = row['canonical_order']
            
            if isinstance(gene_order, str):
                try:
                    gene_order = eval(gene_order)
                except:
                    continue
            
            is_canonical = (gene_order == CANONICAL_ORDER or gene_order == CANONICAL_ORDER[::-1])
            
            if marked_canonical != is_canonical:
                print(f"  Discrepancy in {genome_id}:")
                print(f"    Marked as canonical: {marked_canonical}")
                print(f"    Actually canonical: {is_canonical}")
                print(f"    Gene order: {gene_order}")
    else:
        print("✓ No discrepancies found - analysis appears correct!")

if __name__ == "__main__":
    main()