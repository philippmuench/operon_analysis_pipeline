#!/usr/bin/env python3
"""Test GFF parser on a single genome."""

def test_parse():
    gff_file = "../01_prokka_annotation/output/prokka_results/ENT_AA0002AA_AS.result/ENT_AA0002AA_AS.result.gff"
    
    genes = set()
    gene_count = 0
    total_lines = 0
    
    with open(gff_file, 'r') as f:
        for line in f:
            total_lines += 1
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 9:
                if parts[2] == 'gene':
                    gene_count += 1
                    attrs = parts[8]
                    print(f"Line {total_lines}: {attrs}")
                    for attr in attrs.split(';'):
                        if attr.startswith('gene='):
                            gene_name = attr.replace('gene=', '').strip()
                            genes.add(gene_name)
                            print(f"  Found gene: {gene_name}")
                            break
                    if gene_count >= 5:  # Just show first 5
                        break
    
    print(f"\nTotal genes found: {len(genes)}")
    print(f"Gene names: {sorted(genes)}")

if __name__ == '__main__':
    test_parse()