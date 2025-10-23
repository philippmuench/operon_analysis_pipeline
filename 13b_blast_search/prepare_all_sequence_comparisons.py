#!/usr/bin/env python3
"""
Process all genomes to extract sequences and compare to reference sequences.
Single script to run on SLURM with adequate resources.
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
from typing import Dict, List

class SequenceComparatorAll:
    def __init__(self):
        self.blast_dir = Path("output/blast_results")
        self.assembly_dir = Path("../../Efs_assemblies")
        self.reference_dir = Path("../13a_new_reference_operon_extraction/output")
        self.output_dir = Path("output")
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
            print(f"Error extracting sequence from {assembly_file}: {e}")
            
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
    
    def process_genome(self, genome_id: str, verbose: bool = False) -> Dict:
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
            if verbose:
                print(f"Error processing {genome_id}: {e}")
            
        return genome_data
    
    def process_all_genomes(self):
        """Process all genomes with BLAST results."""
        import time
        from datetime import datetime
        
        # Get all BLAST result files
        all_blast_files = sorted(list(self.blast_dir.glob("*_genes_blast.txt")))
        
        print(f"Found {len(all_blast_files)} genomes with BLAST results")
        print(f"Starting processing at {datetime.now()}")
        print("="*60)
        
        all_genome_data = []
        start_time = time.time()
        genomes_with_matches = 0
        genomes_without_matches = 0
        
        for i, blast_file in enumerate(all_blast_files):
            genome_id = blast_file.stem.replace('_genes_blast', '')
            
            # More frequent logging
            if (i + 1) % 10 == 0:
                elapsed = time.time() - start_time
                rate = (i + 1) / elapsed
                remaining = (len(all_blast_files) - (i + 1)) / rate
                print(f"  [{datetime.now().strftime('%H:%M:%S')}] Processed {i + 1}/{len(all_blast_files)} genomes "
                      f"({(i + 1) / len(all_blast_files) * 100:.1f}%) - "
                      f"Rate: {rate:.1f} genomes/sec - "
                      f"Est. remaining: {remaining/60:.1f} min")
                print(f"    Current genome: {genome_id}")
                print(f"    Genomes with matches: {genomes_with_matches}, without: {genomes_without_matches}")
                sys.stdout.flush()  # Force output to be written
            
            genome_data = self.process_genome(genome_id)
            if genome_data['genes']:
                all_genome_data.append(genome_data)
                genomes_with_matches += 1
            else:
                genomes_without_matches += 1
        
        # Final summary
        total_time = time.time() - start_time
        print("\n" + "="*60)
        print(f"PROCESSING COMPLETE at {datetime.now()}")
        print(f"Total time: {total_time/60:.1f} minutes")
        print(f"Average rate: {len(all_blast_files)/total_time:.2f} genomes/sec")
        print(f"Genomes with matches: {genomes_with_matches}")
        print(f"Genomes without matches: {genomes_without_matches}")
        print("="*60)
        
        print(f"\nSaving results to JSON...")
        sys.stdout.flush()
        
        # Save results
        output_file = self.output_dir / "sequence_variants.json"
        with open(output_file, 'w') as f:
            json.dump(all_genome_data, f, indent=2)
            
        print(f"Results saved to {output_file}")
        
        # Generate summary statistics
        self.generate_summary(all_genome_data)
        
    def generate_summary(self, genome_data: List[Dict]):
        """Generate summary statistics of variants."""
        print("\n" + "="*60)
        print("VARIANT SUMMARY")
        print("="*60)
        
        gene_stats = defaultdict(lambda: {
            'total_snps': 0,
            'total_insertions': 0,
            'total_deletions': 0,
            'avg_identity': [],
            'genomes_analyzed': 0
        })
        
        for genome in genome_data:
            for gene_id, gene_data in genome['genes'].items():
                gene_name = gene_data['gene_name']
                gene_stats[gene_name]['total_snps'] += len(gene_data['snps'])
                gene_stats[gene_name]['total_insertions'] += len(gene_data['insertions'])
                gene_stats[gene_name]['total_deletions'] += len(gene_data['deletions'])
                gene_stats[gene_name]['avg_identity'].append(gene_data['percent_identity'])
                gene_stats[gene_name]['genomes_analyzed'] += 1
        
        # Print summary
        for gene_name, stats in sorted(gene_stats.items()):
            print(f"\n{gene_name}:")
            print(f"  Genomes analyzed: {stats['genomes_analyzed']}")
            print(f"  Total SNPs: {stats['total_snps']}")
            print(f"  Total insertions: {stats['total_insertions']}")
            print(f"  Total deletions: {stats['total_deletions']}")
            if stats['avg_identity']:
                avg_id = np.mean(stats['avg_identity'])
                print(f"  Average identity: {avg_id:.2f}%")
        
        # Save summary to CSV
        summary_df = pd.DataFrame.from_dict(gene_stats, orient='index')
        summary_df['avg_identity'] = summary_df['avg_identity'].apply(
            lambda x: np.mean(x) if x else 0
        )
        summary_df.to_csv(self.output_dir / "variant_summary.csv")
        print(f"\nSummary saved to {self.output_dir / 'variant_summary.csv'}")

def main():
    """Main function."""
    comparator = SequenceComparatorAll()
    comparator.process_all_genomes()
    
if __name__ == "__main__":
    main()
