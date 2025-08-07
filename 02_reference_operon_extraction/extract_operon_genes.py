#!/usr/bin/env python3
"""
Extract operon genes from GenBank file for use as query sequences.
"""

from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import os

def parse_genbank(gb_file):
    """Parse GenBank file and extract gene information."""
    genes = []
    
    # Define the operon genes by location and primary annotation
    # Based on the GenBank file analysis, these are the 7 genes in order
    gene_definitions = [
        {'start': 74, 'end': 1078, 'gene': 'frpC', 'product': 'fructoselysine-6-phosphate deglycase', 'type': 'CDS'},
        {'start': 1072, 'end': 2148, 'gene': 'glpC', 'product': 'glucoselysine-6-phosphate deglycase', 'type': 'CDS'},
        {'start': 2159, 'end': 2992, 'gene': 'ptsD', 'product': 'PTS system mannose/fructose/sorbose family IID component', 'type': 'CDS'},
        {'start': 2973, 'end': 3764, 'gene': 'ptsC', 'product': 'PTS system sorbose-specific IIC component', 'type': 'CDS'},
        {'start': 3780, 'end': 4250, 'gene': 'ptsB', 'product': 'PTS system sorbose subfamily IIB component', 'type': 'CDS'},
        {'start': 4262, 'end': 4726, 'gene': 'ptsA', 'product': 'PTS system fructose IIA component', 'type': 'CDS'},
        {'start': 4753, 'end': 7512, 'gene': 'fruR', 'product': 'sigma-54 dependent transcriptional regulator', 'type': 'CDS'},
        {'start': 7513, 'end': None, 'gene': 'promoter', 'product': 'operon promoter region', 'type': 'promoter'},
        {'start': 7555, 'end': 7562, 'gene': 'pribnow_box', 'product': 'Pribnow box (-10 box)', 'type': 'regulatory'}
    ]
    
    with open(gb_file, 'r') as f:
        for record in SeqIO.parse(f, 'genbank'):
            sequence = str(record.seq)
            
            # Extract each gene based on known coordinates
            for i, gene_def in enumerate(gene_definitions, 1):
                start = gene_def['start'] - 1  # Convert to 0-based
                end = gene_def['end'] if gene_def['end'] else len(sequence)
                
                # All features in this operon are on the complement strand
                gene_seq = sequence[start:end]
                gene_seq_rc = Seq(gene_seq).reverse_complement()
                
                # For CDS features, translate to protein
                if gene_def['type'] == 'CDS':
                    try:
                        protein_seq = str(gene_seq_rc.translate(to_stop=True))
                    except:
                        protein_seq = str(gene_seq_rc.translate())
                else:
                    # For promoter, keep DNA sequence (reverse complement)
                    protein_seq = str(gene_seq_rc)
                
                # Determine locus tag based on feature type
                if gene_def['type'] == 'CDS':
                    cds_count = len([g for g in gene_definitions[:i-1] if g["type"] == "CDS"]) + 1
                    locus_tag = f'operon_CDS_{cds_count}'
                elif gene_def['type'] == 'promoter':
                    locus_tag = 'operon_promoter'
                else:
                    locus_tag = f'operon_{gene_def["gene"]}'
                
                gene_info = {
                    'locus_tag': locus_tag,
                    'gene': gene_def['gene'],
                    'product': gene_def['product'],
                    'start': gene_def['start'],
                    'end': end,  # Use calculated end
                    'strand': -1,  # All on complement strand
                    'translation': protein_seq
                }
                genes.append(gene_info)
    
    return pd.DataFrame(genes)

def main():
    # Input and output files
    gb_file = 'operon.gb'  # Now using local file
    output_dir = 'output'
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Parse GenBank file
    print(f"Parsing {gb_file}...")
    df = parse_genbank(gb_file)
    
    # Save gene table
    output_tsv = os.path.join(output_dir, 'operon_genes.tsv')
    df.to_csv(output_tsv, sep='\t', index=False)
    print(f"Saved gene information to {output_tsv}")
    
    # Save non-coding sequences (promoter, Pribnow box) in nucleotide FASTA format
    output_noncoding_fasta = os.path.join(output_dir, 'operon_noncoding_nt.fasta')
    with open(output_noncoding_fasta, 'w') as f:
        for record in SeqIO.parse(gb_file, 'genbank'):
            sequence = str(record.seq)
            
            for _, gene in df.iterrows():
                if gene['gene'] in ['promoter', 'pribnow_box']:
                    # Extract nucleotide sequence for non-coding features
                    start = gene['start'] - 1  # Convert to 0-based
                    end = gene['end']
                    
                    # Get DNA sequence (reverse complement as all features are on minus strand)
                    gene_seq = sequence[start:end]
                    gene_seq = str(Seq(gene_seq).reverse_complement())
                    
                    # Create descriptive header
                    header = f"{gene['locus_tag']}|{gene['gene']}|{gene['product']}"
                    header = header.replace(' ', '_').replace(',', '')
                    f.write(f">{header}\n{gene_seq}\n")
    
    print(f"Saved non-coding sequences to {output_noncoding_fasta}")
    
    # Save gene sequences in nucleotide FASTA format (for BLASTP after translation)
    output_genes_nt_fasta = os.path.join(output_dir, 'operon_genes_nt.fasta')
    with open(output_genes_nt_fasta, 'w') as f:
        for record in SeqIO.parse(gb_file, 'genbank'):
            sequence = str(record.seq)
            
            for _, gene in df.iterrows():
                if gene['gene'] not in ['promoter', 'pribnow_box']:
                    # Extract nucleotide sequence for CDS features
                    start = gene['start'] - 1  # Convert to 0-based
                    end = gene['end']
                    
                    # Get DNA sequence (reverse complement as all features are on minus strand)
                    gene_seq = sequence[start:end]
                    gene_seq = str(Seq(gene_seq).reverse_complement())
                    
                    # Create descriptive header
                    header = f"{gene['locus_tag']}|{gene['gene']}|{gene['product']}"
                    header = header.replace(' ', '_').replace(',', '')
                    f.write(f">{header}\n{gene_seq}\n")
    
    print(f"Saved gene nucleotide sequences to {output_genes_nt_fasta}")
    
    # Save protein sequences for gene BLAST searches
    output_genes_protein_fasta = os.path.join(output_dir, 'operon_genes_protein.fasta')
    with open(output_genes_protein_fasta, 'w') as f:
        for _, gene in df.iterrows():
            if gene['gene'] not in ['promoter', 'pribnow_box'] and gene['translation']:
                # Create descriptive header
                header = f"{gene['locus_tag']}|{gene['gene']}|{gene['product']}"
                header = header.replace(' ', '_').replace(',', '')
                f.write(f">{header}\n{gene['translation']}\n")
    
    print(f"Saved gene protein sequences to {output_genes_protein_fasta}")
    
    # Display summary
    print(f"\nExtracted {len(df)} genes from operon")
    print("\nGene summary:")
    for _, gene in df.iterrows():
        print(f"  {gene['gene']}: {gene['product']}")

if __name__ == '__main__':
    main()