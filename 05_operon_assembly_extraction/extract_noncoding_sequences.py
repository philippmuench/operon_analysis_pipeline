#!/usr/bin/env python3
"""
Extract non-coding operon sequences (promoter) from genomes that have them.
Uses BLAST results to find and extract promoter sequences.
"""

import os
import sys
import time
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

def extract_noncoding_sequences(prokka_dir, blast_results_dir, output_dir, min_identity=90, min_coverage=50):
    """Extract promoter sequences from genomes."""
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Get list of genomes with complete operons from summary
    summary_file = "../03_blast_search/output/operon_simple_summary.csv"
    if os.path.exists(summary_file):
        summary_df = pd.read_csv(summary_file)
        # Get genomes that have promoters
        promoter_genomes = summary_df[summary_df['has_promoter'] == 1]['genome_id'].tolist()
        print(f"Found {len(promoter_genomes)} genomes with promoters")
    else:
        print("Warning: No summary file found, processing all genomes with BLAST hits")
        promoter_genomes = None
    
    # Non-coding elements to extract
    noncoding_elements = ['promoter']
    
    # Store sequences for each element
    noncoding_sequences = defaultdict(list)
    
    # Process each genome
    blast_files = sorted([f for f in os.listdir(blast_results_dir) if f.endswith('_noncoding_blast.txt')])
    total_genomes = len(blast_files)
    print(f"Found {total_genomes} genomes with *_noncoding_blast.txt in {blast_results_dir}")
    sys.stdout.flush()
    
    start_time = time.time()
    for index, blast_file in enumerate(blast_files, start=1):
        genome_id = blast_file.replace('_noncoding_blast.txt', '')
        if index == 1 or index % 25 == 0 or index == total_genomes:
            elapsed = time.time() - start_time
            rate = index / elapsed if elapsed > 0 else 0
            remaining = (total_genomes - index) / rate if rate > 0 else 0
            print(f"[{index}/{total_genomes}] Processing {genome_id} (elapsed {elapsed/60:.1f}m, ETA {remaining/60:.1f}m)")
            sys.stdout.flush()
        
        # Skip if not a genome with promoter (if we have the info)
        if promoter_genomes and genome_id not in promoter_genomes:
            continue
        
        # Read BLAST results
        blast_path = os.path.join(blast_results_dir, blast_file)
        if os.path.getsize(blast_path) == 0:
            continue
        
        # Parse BLAST hits
        best_hits = {}
        with open(blast_path, 'r') as f:
            for line in f:
                if line.strip():  # Skip empty lines
                    hit = parse_blast_hit(line)
                    element_name = hit['qseqid'].split('|')[1]
                    
                    # Keep only high-quality hits
                    if hit['pident'] >= min_identity and hit['qcovs'] >= min_coverage:
                        # Keep best hit per element
                        if element_name not in best_hits or hit['pident'] > best_hits[element_name]['pident']:
                            best_hits[element_name] = hit
        
        # Skip if no promoter found
        if 'promoter' not in best_hits:
            continue
        
        # Load genome sequence - try different possible paths
        possible_paths = [
            os.path.join(prokka_dir, f"{genome_id}/prokka_output/{genome_id}.fna"),
            os.path.join(prokka_dir, f"{genome_id}/{genome_id}.fna"),
            os.path.join(prokka_dir, f"{genome_id}.fna")
        ]
        
        genome_file = None
        for path in possible_paths:
            if os.path.exists(path):
                genome_file = path
                break
        
        if not genome_file:
            print(f"Warning: Genome file not found for {genome_id}")
            continue
        else:
            print(f"  Using genome FASTA: {genome_file}")
            sys.stdout.flush()
        
        # Read genome sequences
        genome_seqs = {}
        for record in SeqIO.parse(genome_file, "fasta"):
            genome_seqs[record.id] = record.seq
        
        # Extract sequences
        for element in noncoding_elements:
            if element in best_hits:
                hit = best_hits[element]
                
                # Find sequence with improved ID matching
                seq_found = False
                
                # Extract core ID from BLAST hit (remove gnl|X| prefix and get base name)
                blast_id = hit['sseqid']
                if 'gnl|X|' in blast_id:
                    blast_core = blast_id.split('gnl|X|')[1]
                    # Get the base name without suffix numbers
                    if '_' in blast_core:
                        blast_base = '_'.join(blast_core.split('_')[:-1])
                    else:
                        blast_base = blast_core.rstrip('0123456789')
                else:
                    blast_base = blast_id
                
                for seq_id, seq in genome_seqs.items():
                    # Try multiple matching strategies
                    match_found = False
                    
                    # Strategy 1: Direct substring match (original)
                    if hit['sseqid'] in seq_id or seq_id in hit['sseqid']:
                        match_found = True
                    
                    # Strategy 2: Core ID matching
                    elif blast_core in seq_id or seq_id in blast_core:
                        match_found = True
                    
                    # Strategy 3: Base name matching (handles BEOEMDHG_1 vs BEOEMDHG_15)
                    elif blast_base in seq_id:
                        match_found = True
                    
                    # Strategy 4: Just use the first/longest sequence if no match
                    # (many assemblies have single contigs)
                    
                    if match_found:
                        print(f"  Found sequence match: {hit['sseqid']} -> {seq_id}")
                        
                        # Extract sequence with bounds checking
                        start = min(hit['sstart'], hit['send']) - 1  # Convert to 0-based
                        end = max(hit['sstart'], hit['send'])
                        
                        # Check if coordinates are within sequence bounds
                        if end <= len(seq) and start >= 0:
                            extracted_seq = seq[start:end]
                            
                            # Handle reverse complement if needed
                            if hit['sstart'] > hit['send']:
                                extracted_seq = extracted_seq.reverse_complement()
                            
                            # Only add if we got a reasonable sequence
                            if len(extracted_seq) > 0:
                                # Create sequence record
                                seq_record = f">{genome_id}.result_{element}\n{extracted_seq}\n"
                                noncoding_sequences[element].append(seq_record)
                                seq_found = True
                                print(f"    Successfully extracted {len(extracted_seq)} bp")
                                break
                            else:
                                print(f"    Warning: Extracted sequence is empty")
                        else:
                            print(f"    Warning: Coordinates out of bounds - start:{start+1}, end:{end}, seq_len:{len(seq)}")
                            # Try to find the sequence in other contigs
                            continue
                
                # Fallback: use the longest sequence if no ID match found
                if not seq_found and genome_seqs:
                    longest_seq_id = max(genome_seqs.keys(), key=lambda x: len(genome_seqs[x]))
                    longest_seq = genome_seqs[longest_seq_id]
                    
                    print(f"  No ID match found for {hit['sseqid']}, using longest sequence: {longest_seq_id}")
                    
                    # Extract sequence
                    start = min(hit['sstart'], hit['send']) - 1  # Convert to 0-based
                    end = max(hit['sstart'], hit['send'])
                    
                    if end <= len(longest_seq) and start >= 0:
                        extracted_seq = longest_seq[start:end]
                        
                        # Handle reverse complement if needed
                        if hit['sstart'] > hit['send']:
                            extracted_seq = extracted_seq.reverse_complement()
                        
                        # Only add if we got a reasonable sequence
                        if len(extracted_seq) > 0:
                            # Create sequence record
                            seq_record = f">{genome_id}.result_{element}\n{extracted_seq}\n"
                            noncoding_sequences[element].append(seq_record)
                            seq_found = True
                            print(f"    Successfully extracted {len(extracted_seq)} bp from longest sequence")
                        else:
                            print(f"    Warning: Extracted sequence from longest sequence is empty")
                    else:
                        print(f"    Warning: Coordinates still out of bounds for longest sequence - start:{start+1}, end:{end}, seq_len:{len(longest_seq)}")
                
                if not seq_found:
                    print(f"Warning: Could not find sequence for {element} in {genome_id}")
        print(f"  Extracted promoter for {genome_id}: {1 if 'promoter' in best_hits and seq_found else 0}")
        sys.stdout.flush()
    
    # Write sequences to files
    for element in noncoding_elements:
        if noncoding_sequences[element]:
            output_file = os.path.join(output_dir, f"{element}.fasta")
            with open(output_file, 'w') as f:
                f.writelines(noncoding_sequences[element])
            print(f"Wrote {len(noncoding_sequences[element])} {element} sequences to {output_file}")
        else:
            print(f"Warning: No {element} sequences found")
    
    return {element: len(noncoding_sequences[element]) for element in noncoding_elements}

def main():
    parser = argparse.ArgumentParser(description="Extract non-coding operon sequences")
    parser.add_argument("--prokka-dir", default="../01_prokka_annotation/output/prokka_results", 
                        help="Directory containing Prokka results")
    parser.add_argument("--blast-dir", default="../03_blast_search/output/blast_results", 
                        help="Directory containing BLAST results")
    parser.add_argument("--output-dir", default="output/noncoding_sequences", 
                        help="Output directory for sequences")
    parser.add_argument("--min-identity", type=float, default=90, 
                        help="Minimum sequence identity")
    parser.add_argument("--min-coverage", type=float, default=50, 
                        help="Minimum query coverage")
    
    args = parser.parse_args()
    
    # Extract sequences
    counts = extract_noncoding_sequences(
        args.prokka_dir, 
        args.blast_dir, 
        args.output_dir,
        args.min_identity,
        args.min_coverage
    )
    
    print("\nExtraction summary:")
    for element, count in counts.items():
        print(f"  {element}: {count} sequences")

if __name__ == "__main__":
    main()