#!/usr/bin/env python3
"""
Extract operon gene sequences from genomes that have them.
Uses BLAST results to find and extract sequences.
"""

import os
import sys
import time
import gzip
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import argparse

def parse_blast_hit(blast_line):
    """Parse a BLAST hit line."""
    fields = blast_line.strip().split('\t')
    return {
        'qseqid': fields[0],
        'sseqid': fields[1],
        'pident': float(fields[2]),
        'length': int(fields[3]),
        'sstart': int(fields[8]),
        'send': int(fields[9]),
        'qcovs': float(fields[12]) if len(fields) > 12 else 100.0
    }

def _find_fna_file(prokka_dir: str, genome_id: str):
    """Return path to the .fna file for a given genome_id, handling .result suffixes.

    Tries common layouts produced by Prokka:
    - {prokka_dir}/{genome_id}.result/{genome_id}.result.fna
    - {prokka_dir}/{genome_id}.result/{genome_id}.fna
    - {prokka_dir}/{genome_id}/{genome_id}.fna
    - {prokka_dir}/{genome_id}.fna
    - Any single *.fna inside the matched directory
    """
    # Candidate directories
    candidate_dirs = [
        os.path.join(prokka_dir, f"{genome_id}.result"),
        os.path.join(prokka_dir, genome_id),
    ]

    # If there is exactly one dir that starts with genome_id and ends with .result, use it
    try:
        for entry in os.listdir(prokka_dir):
            full = os.path.join(prokka_dir, entry)
            if os.path.isdir(full) and entry.startswith(genome_id) and entry.endswith('.result'):
                if full not in candidate_dirs:
                    candidate_dirs.insert(0, full)
    except FileNotFoundError:
        pass

    # Try each directory for expected files
    for d in candidate_dirs:
        if not os.path.isdir(d):
            continue
        dir_basename = os.path.basename(d)
        fna_candidates = [
            os.path.join(d, f"{dir_basename}.fna"),
            os.path.join(d, f"{genome_id}.fna"),
        ]
        for cand in fna_candidates:
            if os.path.exists(cand):
                return cand
        # Fallback: any single .fna in directory
        try:
            fna_files = [os.path.join(d, f) for f in os.listdir(d) if f.endswith('.fna')]
            if len(fna_files) == 1:
                return fna_files[0]
            elif len(fna_files) > 1:
                # Prefer one that contains genome_id
                for f in fna_files:
                    if genome_id in os.path.basename(f):
                        return f
                # Else return the first
                return fna_files[0]
        except FileNotFoundError:
            pass

    # As a last resort, look for a top-level file
    top_level = os.path.join(prokka_dir, f"{genome_id}.fna")
    if os.path.exists(top_level):
        return top_level
    return None


def _find_assembly_fasta(assemblies_dir: str, genome_id: str):
    """Return path to the raw assembly FASTA for a given genome_id.

    Tries exact filenames only, matching the provided genome_id without
    altering it (genome_id already includes potential suffixes like
    '.result'):
    - {assemblies_dir}/{genome_id}.fasta.gz
    - {assemblies_dir}/{genome_id}.fasta
    """
    if not assemblies_dir:
        return None
    candidates = [
        os.path.join(assemblies_dir, f"{genome_id}.fasta.gz"),
        os.path.join(assemblies_dir, f"{genome_id}.fasta"),
    ]
    # Backward-compatible fallbacks if callers passed a genome_id without
    # the '.result' suffix but assemblies are stored with it
    if not genome_id.endswith('.result'):
        candidates.extend([
            os.path.join(assemblies_dir, f"{genome_id}.result.fasta.gz"),
            os.path.join(assemblies_dir, f"{genome_id}.result.fasta"),
        ])
    for cand in candidates:
        if os.path.exists(cand):
            return cand
    return None


def _read_genome_sequences(genome_file: str):
    """Read sequences from a FASTA file (supports .gz) and return id->record mapping."""
    sequences = {}
    handle = gzip.open(genome_file, 'rt') if genome_file.endswith('.gz') else open(genome_file, 'r')
    try:
        for record in SeqIO.parse(handle, "fasta"):
            sequences[record.id] = record
    finally:
        handle.close()
    return sequences


def extract_sequences(prokka_dir, blast_results_dir, output_dir, min_identity=90, min_coverage=80, source="prokka", assemblies_dir=None):
    """Extract operon gene sequences from genomes."""
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Get list of genomes with complete operons from summary
    summary_file = "../03_blast_search/output/operon_simple_summary.csv"
    if os.path.exists(summary_file):
        summary_df = pd.read_csv(summary_file)
        complete_genomes = summary_df[summary_df['is_complete'] == 1]['genome_id'].tolist()
        print(f"Found {len(complete_genomes)} genomes with complete operons")
    else:
        print("Warning: No summary file found, processing all genomes with BLAST hits")
        complete_genomes = None
    
    # Operon genes to extract
    operon_genes = ['frpC', 'glpC', 'ptsD', 'ptsC', 'ptsB', 'ptsA', 'fruR']
    
    # Store sequences for each gene
    gene_sequences = defaultdict(list)
    
    # Process each genome
    blast_files = sorted([f for f in os.listdir(blast_results_dir) if f.endswith('_genes_blast.txt')])
    total_genomes = len(blast_files)
    print(f"Found {total_genomes} genomes with *_genes_blast.txt in {blast_results_dir}")
    sys.stdout.flush()
    
    start_time = time.time()
    for index, blast_file in enumerate(blast_files, start=1):
        genome_id = blast_file.replace('_genes_blast.txt', '')
        if index == 1 or index % 25 == 0 or index == total_genomes:
            elapsed = time.time() - start_time
            rate = index / elapsed if elapsed > 0 else 0
            remaining = (total_genomes - index) / rate if rate > 0 else 0
            print(f"[{index}/{total_genomes}] Processing {genome_id} (elapsed {elapsed/60:.1f}m, ETA {remaining/60:.1f}m)")
            sys.stdout.flush()
        
        # Skip if not a complete genome (if we have the info)
        if complete_genomes and genome_id not in complete_genomes:
            continue
        
        # Read BLAST results
        blast_path = os.path.join(blast_results_dir, blast_file)
        if os.path.getsize(blast_path) == 0:
            continue
        
        # Parse BLAST hits
        best_hits = {}
        with open(blast_path, 'r') as f:
            for line in f:
                hit = parse_blast_hit(line)
                gene_name = hit['qseqid'].split('|')[1]
                
                # Keep only high-quality hits
                if hit['pident'] >= min_identity and hit['qcovs'] >= min_coverage:
                    # Keep best hit per gene
                    if gene_name not in best_hits or hit['pident'] > best_hits[gene_name]['pident']:
                        best_hits[gene_name] = hit
        print(f"  {genome_id}: candidate genes passing filters: {len(best_hits)}/{len(operon_genes)}")
        sys.stdout.flush()
        
        # Check if we have all operon genes
        if not all(gene in best_hits for gene in operon_genes):
            continue
        
        # Load genome sequence from requested source (no cross-source fallback)
        if source == "prokka":
            genome_file = _find_fna_file(prokka_dir, genome_id)
        elif source == "assemblies":
            genome_file = _find_assembly_fasta(assemblies_dir, genome_id)
        else:
            genome_file = None

        if not genome_file:
            print(f"Warning: Genome file not found for {genome_id} (source={source})")
            continue
        else:
            print(f"  Using genome FASTA: {genome_file}")
            sys.stdout.flush()
        
        # Extract sequences
        sequences = {}
        records_by_id = _read_genome_sequences(genome_file)
        for gene, hit in best_hits.items():
            if gene not in operon_genes:
                continue
            record = records_by_id.get(hit['sseqid'])
            if record is None:
                continue
            # Extract sequence based on coordinates
            start = min(hit['sstart'], hit['send']) - 1
            end = max(hit['sstart'], hit['send'])

            if hit['sstart'] > hit['send']:
                seq = record.seq[start:end].reverse_complement()
            else:
                seq = record.seq[start:end]

            sequences[gene] = seq
        
        # Add sequences to collection
        added_this_genome = 0
        for gene in operon_genes:
            if gene in sequences:
                gene_sequences[gene].append({
                    'genome_id': genome_id,
                    'sequence': str(sequences[gene])
                })
                added_this_genome += 1
        print(f"  Extracted {added_this_genome}/{len(operon_genes)} operon genes for {genome_id}")
        sys.stdout.flush()
    
    # Write sequences to files
    # Ensure output directory still exists before writing (extra safety in case it was removed mid-run)
    os.makedirs(output_dir, exist_ok=True)
    for gene, seqs in gene_sequences.items():
        output_file = os.path.join(output_dir, f"{gene}.fasta")
        # Double-check parent directory exists (handles unusual path configurations)
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        with open(output_file, 'w') as f:
            for seq_info in seqs:
                f.write(f">{seq_info['genome_id']}\n{seq_info['sequence']}\n")
        print(f"Wrote {len(seqs)} sequences for {gene} to {output_file}")
        sys.stdout.flush()
    
    # Create summary
    summary = {
        'total_genomes': len(complete_genomes) if complete_genomes else len(blast_files),
        'genomes_processed': len(set(s['genome_id'] for seqs in gene_sequences.values() for s in seqs)),
        'genes_extracted': list(gene_sequences.keys()),
        'sequences_per_gene': {gene: len(seqs) for gene, seqs in gene_sequences.items()}
    }
    
    return summary

def main():
    parser = argparse.ArgumentParser(description="Extract operon sequences from genomes")
    parser.add_argument('--prokka_dir', default='../01_prokka_annotation/output/prokka_results',
                        help='Directory with Prokka outputs')
    parser.add_argument('--blast_dir', default='../03_blast_search/output/blast_results',
                        help='Directory with BLAST results')
    parser.add_argument('--output_dir', default='output/operon_sequences',
                        help='Output directory for sequences')
    parser.add_argument('--min_identity', type=float, default=90,
                        help='Minimum identity threshold')
    parser.add_argument('--min_coverage', type=float, default=80,
                        help='Minimum coverage threshold')
    parser.add_argument('--source', choices=['prokka', 'assemblies'], default='prokka',
                        help="Choose genome source: 'prokka' for Prokka .fna, 'assemblies' for raw assemblies")
    parser.add_argument('--assemblies_dir', default='../Efs_assemblies',
                        help='Directory containing raw assemblies when --source=assemblies')
    
    args = parser.parse_args()
    
    print("Extracting operon sequences...")
    summary = extract_sequences(
        args.prokka_dir,
        args.blast_dir,
        args.output_dir,
        args.min_identity,
        args.min_coverage,
        args.source,
        args.assemblies_dir,
    )
    
    print(f"\nExtraction complete:")
    print(f"Total genomes: {summary['total_genomes']}")
    print(f"Genomes with complete operons: {summary['genomes_processed']}")
    print(f"Sequences per gene:")
    for gene, count in summary['sequences_per_gene'].items():
        print(f"  {gene}: {count}")

if __name__ == '__main__':
    main()