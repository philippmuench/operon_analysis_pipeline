#!/usr/bin/env python3
"""
Extract nucleotide sequences for core genes from Prokka output.
Creates FASTA files for each core gene across all genomes.
"""

import os
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from multiprocessing import Pool, cpu_count
import argparse

def extract_gene_from_genome(args):
    """Extract a specific gene from a genome's Prokka output."""
    gene_name, genome_dir = args
    
    genome_id = os.path.basename(genome_dir.rstrip('/'))
    
    # Files to check
    gff_file = os.path.join(genome_dir, f"{genome_id}.gff")
    fna_file = os.path.join(genome_dir, f"{genome_id}.fna")
    
    if not os.path.exists(gff_file) or not os.path.exists(fna_file):
        return None
    
    # Parse genome sequence
    genome_seqs = {}
    try:
        for record in SeqIO.parse(fna_file, "fasta"):
            genome_seqs[record.id] = record
    except:
        return None
    
    # Find gene coordinates in GFF
    gene_coords = None
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 9 and parts[2] == 'gene':
                # Check if this is our gene
                attrs = parts[8]
                for attr in attrs.split(';'):
                    if attr.startswith('gene=') and attr.replace('gene=', '').strip() == gene_name:
                        gene_coords = {
                            'contig': parts[0],
                            'start': int(parts[3]),
                            'end': int(parts[4]),
                            'strand': parts[6]
                        }
                        break
                if gene_coords:
                    break
    
    if not gene_coords:
        return None
    
    # Extract sequence
    contig_id = gene_coords['contig']
    if contig_id not in genome_seqs:
        return None
    
    contig_seq = genome_seqs[contig_id].seq
    start = gene_coords['start'] - 1  # Convert to 0-based
    end = gene_coords['end']
    
    gene_seq = contig_seq[start:end]
    
    # Reverse complement if on negative strand
    if gene_coords['strand'] == '-':
        gene_seq = gene_seq.reverse_complement()
    
    # Create SeqRecord
    seq_record = SeqRecord(
        gene_seq,
        id=f"{genome_id}_{gene_name}",
        description=f"{gene_name} from {genome_id} [{gene_coords['contig']}:{start+1}-{end}({gene_coords['strand']})]"
    )
    
    return gene_name, genome_id, seq_record

def main():
    parser = argparse.ArgumentParser(description="Extract core gene sequences from Prokka output")
    parser.add_argument("--core-genes", default="output/core_genes_95pct.txt",
                        help="File containing list of core genes")
    parser.add_argument("--prokka-dir", default="../01_prokka_annotation/output/prokka_results",
                        help="Directory containing Prokka results")
    parser.add_argument("--output-dir", default="output/core_gene_sequences",
                        help="Output directory for gene sequences")
    parser.add_argument("--threads", type=int, default=min(cpu_count()-1, 20),
                        help="Number of threads to use")
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load core genes
    if not os.path.exists(args.core_genes):
        print(f"Error: Core genes file not found: {args.core_genes}")
        print("Run identify_core_genes.py first!")
        return 1
    
    with open(args.core_genes, 'r') as f:
        core_genes = [line.strip() for line in f if line.strip()]
    
    print(f"Loaded {len(core_genes)} core genes")
    
    # Get genome directories
    genome_dirs = glob.glob(os.path.join(args.prokka_dir, '*/'))
    genome_dirs = [d.rstrip('/') for d in genome_dirs]
    
    print(f"Found {len(genome_dirs)} genomes")
    
    if not genome_dirs:
        print(f"Error: No genomes found in {args.prokka_dir}")
        return 1
    
    # Process each core gene
    total_genes = len(core_genes)
    gene_sequences = {}
    
    for i, gene_name in enumerate(core_genes, 1):
        print(f"Processing gene {i}/{total_genes}: {gene_name}")
        
        # Prepare arguments for parallel processing
        task_args = [(gene_name, genome_dir) for genome_dir in genome_dirs]
        
        # Extract sequences in parallel
        with Pool(args.threads) as pool:
            results = pool.map(extract_gene_from_genome, task_args)
        
        # Collect valid results
        gene_seqs = []
        for result in results:
            if result:
                gene, genome_id, seq_record = result
                gene_seqs.append(seq_record)
        
        if gene_seqs:
            gene_sequences[gene_name] = gene_seqs
            
            # Save sequences for this gene
            output_file = os.path.join(args.output_dir, f"{gene_name}.fasta")
            SeqIO.write(gene_seqs, output_file, "fasta")
            
            print(f"  → Extracted {len(gene_seqs)} sequences for {gene_name}")
        else:
            print(f"  → No sequences found for {gene_name}")
    
    # Summary
    print(f"\n{'='*60}")
    print(f"Core gene sequence extraction complete!")
    print(f"{'='*60}")
    print(f"Total core genes: {len(core_genes)}")
    print(f"Genes with sequences: {len(gene_sequences)}")
    print(f"Output directory: {args.output_dir}")
    
    # Create summary file
    summary_file = os.path.join(args.output_dir, "extraction_summary.txt")
    with open(summary_file, 'w') as f:
        f.write("Core Gene Sequence Extraction Summary\n")
        f.write("=====================================\n\n")
        f.write(f"Total core genes: {len(core_genes)}\n")
        f.write(f"Genes with sequences: {len(gene_sequences)}\n")
        f.write(f"Genomes processed: {len(genome_dirs)}\n\n")
        f.write("Genes successfully extracted:\n")
        for gene_name in sorted(gene_sequences.keys()):
            f.write(f"  {gene_name}: {len(gene_sequences[gene_name])} sequences\n")
    
    print(f"Summary saved to: {summary_file}")
    return 0

if __name__ == "__main__":
    exit(main())