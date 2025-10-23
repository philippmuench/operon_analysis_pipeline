#!/usr/bin/env python3
"""
Batch processing script for sequence comparison.
Processes a subset of genomes based on batch ID.
"""

import os
import sys
import gzip
import json
import pandas as pd
import numpy as np
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import argparse
from typing import Dict, List

class SequenceComparatorBatch:
    def __init__(self, batch_id: int, batch_size: int, output_dir: str):
        self.batch_id = batch_id
        self.batch_size = batch_size
        self.blast_dir = Path("output/blast_results")
        self.assembly_dir = Path("../../Efs_assemblies")
        self.reference_dir = Path("../13a_new_reference_operon_extraction/output")
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Load reference sequences
        self.ref_sequences = {}
        self.load_reference_sequences()
        
    def load_reference_sequences(self):
        """Load reference sequences from FASTA files."""
        print("Loading reference sequences...")
        
        # Load nucleotide sequences
        nt_file = self.reference_dir / "operon_genes_nt.fasta"
        if nt_file.exists():
            for record in SeqIO.parse(nt_file, "fasta"):
                self.ref_sequences[record.id] = str(record.seq)
                
        # Load non-coding sequences
        nc_file = self.reference_dir / "operon_regulatory_nt.fasta"
        if nc_file.exists():
            for record in SeqIO.parse(nc_file, "fasta"):
                self.ref_sequences[record.id] = str(record.seq)
                
        print(f"Loaded {len(self.ref_sequences)} reference sequences")
        
    def extract_sequence_from_assembly(self, assembly_file: Path, contig: str, 
                                      start: int, end: int, strand: str) -> str:
        """Extract sequence from assembly file at given coordinates."""
        sequences = {}
        
        try:
            # Handle compressed files
            if assembly_file.suffix == '.gz':
                with gzip.open(assembly_file, 'rt') as f:
                    for record in SeqIO.parse(f, "fasta"):
                        sequences[record.id] = str(record.seq)
            else:
                for record in SeqIO.parse(assembly_file, "fasta"):
                    sequences[record.id] = str(record.seq)
            
            # Find the contig
            for seq_id, seq in sequences.items():
                if contig in seq_id:
                    # Extract sequence (BLAST uses 1-based coordinates)
                    extracted = seq[start-1:end]
                    
                    # Reverse complement if on reverse strand
                    if strand == '-' or start > end:
                        extracted = str(Seq(extracted).reverse_complement())
                        
                    return extracted
        except Exception as e:
            print(f"Error extracting sequence: {e}")
            
        return ""
    
    def align_sequences(self, ref_seq: str, query_seq: str) -> Dict:
        """Simple alignment to find variants between reference and query."""
        variants = {
            'snps': [],
            'insertions': [],
            'deletions': [],
            'identity_per_position': []
        }
        
        if not ref_seq or not query_seq:
            variants['percent_identity'] = 0
            variants['conservation_windows'] = []
            return variants
        
        # Simple pairwise alignment
        min_len = min(len(ref_seq), len(query_seq))
        
        # Calculate identity per position
        for i in range(min_len):
            if i < len(ref_seq) and i < len(query_seq):
                if ref_seq[i] == query_seq[i]:
                    variants['identity_per_position'].append(1)
                else:
                    variants['identity_per_position'].append(0)
                    # Record SNP
                    variants['snps'].append({
                        'position': i + 1,
                        'ref': ref_seq[i],
                        'alt': query_seq[i]
                    })
        
        # Handle length differences
        if len(ref_seq) > len(query_seq):
            # Deletion in query
            variants['deletions'].append({
                'position': min_len + 1,
                'length': len(ref_seq) - len(query_seq),
                'sequence': ref_seq[min_len:]
            })
        elif len(query_seq) > len(ref_seq):
            # Insertion in query
            variants['insertions'].append({
                'position': min_len + 1,
                'length': len(query_seq) - len(ref_seq),
                'sequence': query_seq[min_len:]
            })
            
        # Calculate overall identity
        if min_len > 0:
            variants['percent_identity'] = (sum(variants['identity_per_position']) / min_len) * 100
        else:
            variants['percent_identity'] = 0
            
        # Calculate conservation windows (10bp windows)
        window_size = 10
        conservation_windows = []
        for i in range(0, len(variants['identity_per_position']), window_size):
            window = variants['identity_per_position'][i:i+window_size]
            if window:
                conservation_windows.append(sum(window) / len(window))
                
        variants['conservation_windows'] = conservation_windows
        
        return variants
    
    def process_genome(self, genome_id: str) -> Dict:
        """Process a single genome to extract sequences and compare to reference."""
        genome_data = {
            'genome_id': genome_id,
            'genes': {}
        }
        
        # Find assembly file
        assembly_file = self.assembly_dir / f"{genome_id}.result.fasta.gz"
        if not assembly_file.exists():
            assembly_file = self.assembly_dir / f"{genome_id}.fasta.gz"
            if not assembly_file.exists():
                return genome_data
        
        # Load BLAST results
        blast_file = self.blast_dir / f"{genome_id}_genes_blast.txt"
        if not blast_file.exists():
            return genome_data
            
        try:
            # Parse BLAST results
            blast_df = pd.read_csv(blast_file, sep='\t', header=None,
                                   names=['qseqid', 'sseqid', 'pident', 'length', 
                                         'mismatch', 'gapopen', 'qstart', 'qend',
                                         'sstart', 'send', 'evalue', 'bitscore', 'qcovs'])
            
            # Process each gene
            for gene_id in self.ref_sequences.keys():
                # Get best hit for this gene
                gene_hits = blast_df[blast_df['qseqid'] == gene_id]
                if gene_hits.empty:
                    continue
                    
                best_hit = gene_hits.iloc[0]  # Already sorted by bitscore
                
                # Extract sequence from assembly
                strand = '+' if best_hit['sstart'] < best_hit['send'] else '-'
                query_seq = self.extract_sequence_from_assembly(
                    assembly_file,
                    best_hit['sseqid'],
                    min(best_hit['sstart'], best_hit['send']),
                    max(best_hit['sstart'], best_hit['send']),
                    strand
                )
                
                if not query_seq:
                    continue
                    
                # Get reference sequence
                ref_seq = self.ref_sequences[gene_id]
                
                # Compare sequences
                variants = self.align_sequences(ref_seq, query_seq)
                
                # Add metadata
                variants['gene_id'] = gene_id
                variants['gene_name'] = gene_id.split('|')[1] if '|' in gene_id else gene_id
                variants['contig'] = best_hit['sseqid']
                variants['start'] = int(best_hit['sstart'])
                variants['end'] = int(best_hit['send'])
                variants['strand'] = strand
                variants['blast_identity'] = float(best_hit['pident'])
                variants['coverage'] = float(best_hit['qcovs'])
                
                genome_data['genes'][gene_id] = variants
                
        except Exception as e:
            print(f"Error processing {genome_id}: {e}")
            
        return genome_data
    
    def process_batch(self):
        """Process a batch of genomes based on batch ID."""
        # Get all BLAST result files
        all_blast_files = sorted(list(self.blast_dir.glob("*_genes_blast.txt")))
        
        # Calculate batch range
        start_idx = (self.batch_id - 1) * self.batch_size
        end_idx = min(start_idx + self.batch_size, len(all_blast_files))
        
        if start_idx >= len(all_blast_files):
            print(f"Batch {self.batch_id} is beyond the range of available files")
            return
        
        batch_files = all_blast_files[start_idx:end_idx]
        
        print(f"Batch {self.batch_id}: Processing genomes {start_idx + 1} to {end_idx}")
        print(f"Processing {len(batch_files)} genomes")
        
        all_genome_data = []
        
        for i, blast_file in enumerate(batch_files):
            genome_id = blast_file.stem.replace('_genes_blast', '')
            
            if (i + 1) % 10 == 0:
                print(f"  Processed {i + 1}/{len(batch_files)} genomes in batch {self.batch_id}")
            
            genome_data = self.process_genome(genome_id)
            if genome_data['genes']:
                all_genome_data.append(genome_data)
        
        # Save batch results
        output_file = self.output_dir / f"variants_batch_{self.batch_id:04d}.json"
        with open(output_file, 'w') as f:
            json.dump(all_genome_data, f, indent=2)
            
        print(f"Batch {self.batch_id}: Saved {len(all_genome_data)} genomes to {output_file}")
        
        # Save batch summary
        summary = {
            'batch_id': self.batch_id,
            'genomes_processed': len(all_genome_data),
            'genome_ids': [g['genome_id'] for g in all_genome_data]
        }
        
        summary_file = self.output_dir / f"summary_batch_{self.batch_id:04d}.json"
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)

def main():
    parser = argparse.ArgumentParser(description='Batch process sequence comparisons')
    parser.add_argument('--batch-id', type=int, required=True, help='Batch ID (1-based)')
    parser.add_argument('--batch-size', type=int, default=100, help='Number of genomes per batch')
    parser.add_argument('--output-dir', default='output/sequence_variants_parts', 
                       help='Output directory for batch results')
    
    args = parser.parse_args()
    
    comparator = SequenceComparatorBatch(
        batch_id=args.batch_id,
        batch_size=args.batch_size,
        output_dir=args.output_dir
    )
    
    comparator.process_batch()

if __name__ == "__main__":
    main()
