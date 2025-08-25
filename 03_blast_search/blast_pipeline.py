#!/usr/bin/env python3
"""
Unified BLAST pipeline for operon analysis.
Combines BLAST search, result processing, and statistics generation.
"""

import os
import sys
import glob
import argparse
import pandas as pd
import subprocess
from pathlib import Path
from typing import List, Dict, Optional
import logging
from datetime import datetime
import numpy as np

class BlastPipeline:
    def __init__(self, 
                 prokka_dir: str = "../01_prokka_annotation/output/prokka_results",
                 assembly_dir: str = "../../Efs_assemblies",
                 reference_dir: str = "../02_reference_operon_extraction/output",
                 output_dir: str = "output",
                 batch_id: Optional[int] = None,
                 batch_size: int = 100,
                 num_threads: int = 60):
        
        self.prokka_dir = Path(prokka_dir)
        self.assembly_dir = Path(assembly_dir)
        self.reference_dir = Path(reference_dir)
        self.output_dir = Path(output_dir)
        self.batch_id = batch_id
        self.batch_size = batch_size
        self.num_threads = num_threads
        
        # Output directories for different BLAST modes
        self.blast_dirs = {
            'coding_protein': self.output_dir / 'blast_results',
            'coding_nt': self.output_dir / 'blast_results_nt',
            'prokka_variants': self.output_dir / 'blast_results_prokka_variants',
            'noncoding': self.output_dir / 'blast_results'
        }
        
        # Reference query files
        self.query_files = {
            'protein': self.reference_dir / 'operon_genes_protein.fasta',
            'nt': self.reference_dir / 'operon_genes_nt.fasta',
            'noncoding': self.reference_dir / 'operon_noncoding_nt.fasta'
        }
        
        # Setup logging
        self.setup_logging()
        
    def setup_logging(self):
        """Setup logging configuration."""
        log_file = self.output_dir / f'blast_pipeline_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)
        
    def create_directories(self):
        """Create all necessary output directories."""
        for dir_path in self.blast_dirs.values():
            dir_path.mkdir(parents=True, exist_ok=True)
        
        # Create reference DB directory
        ref_db_dir = self.output_dir / 'reference_db'
        ref_db_dir.mkdir(parents=True, exist_ok=True)
        
    def get_batch_genomes(self, use_assemblies: bool = True) -> List[Path]:
        """Get list of genomes for the current batch."""
        if use_assemblies:
            # Get assembly files
            all_files = sorted(self.assembly_dir.glob('*.fasta.gz'))
            self.logger.info(f"Found {len(all_files)} assembly files")
        else:
            # Get prokka directories for variants
            all_files = sorted([d for d in self.prokka_dir.iterdir() if d.is_dir()])
            self.logger.info(f"Found {len(all_files)} prokka directories")
        
        if self.batch_id is not None:
            start_idx = (self.batch_id - 1) * self.batch_size
            end_idx = self.batch_id * self.batch_size
            batch_files = all_files[start_idx:end_idx]
            self.logger.info(f"Processing batch {self.batch_id}: {len(batch_files)} items")
            return batch_files
        else:
            self.logger.info(f"Processing all {len(all_files)} items")
            return all_files
            
    def detect_prefix(self, genome_dir: Path) -> Optional[str]:
        """Detect the file prefix for a genome directory."""
        for suffix in ['.fna', '.faa', '.gff']:
            files = list(genome_dir.glob(f'*{suffix}'))
            if files:
                return files[0].stem
        return None
        
    def run_blast_search(self, modes: List[str] = None):
        """Run BLAST searches for specified modes."""
        if modes is None:
            modes = ['coding_protein', 'coding_nt', 'prokka_variants', 'noncoding']
            
        self.logger.info(f"Running BLAST searches for modes: {modes}")
        
        # Handle prokka_variants separately
        if 'prokka_variants' in modes:
            self.run_prokka_variants_search()
            modes = [m for m in modes if m != 'prokka_variants']
        
        # Process assembly-based searches
        if modes:
            batch_genomes = self.get_batch_genomes(use_assemblies=True)
            
            for assembly_file in batch_genomes:
                if not assembly_file.exists():
                    continue
                    
                # Get genome name from assembly file
                genome_name = assembly_file.stem.replace('.result.fasta', '')
                self.logger.info(f"Processing genome: {genome_name}")
                
                # Run each mode against assembly
                if 'coding_protein' in modes:
                    self.run_tblastn(genome_name, assembly_file)
                    
                if 'coding_nt' in modes:
                    self.run_blastn_coding(genome_name, assembly_file)
                    
                if 'noncoding' in modes:
                    self.run_blastn_noncoding(genome_name, assembly_file)
                
    def prepare_blast_db(self):
        """Prepare BLAST database for prokka_variants mode."""
        ref_db_dir = self.output_dir / 'reference_db'
        db_fasta = ref_db_dir / 'operon_genes_nt.fasta'
        
        # Copy reference file
        if self.query_files['nt'].exists():
            subprocess.run(['cp', '-f', str(self.query_files['nt']), str(db_fasta)])
            
            # Make BLAST database if it doesn't exist
            db_index = ref_db_dir / 'operon_genes_nt.nsq'
            if not db_index.exists():
                cmd = [
                    'makeblastdb',
                    '-in', str(db_fasta),
                    '-dbtype', 'nucl',
                    '-out', str(ref_db_dir / 'operon_genes_nt')
                ]
                subprocess.run(cmd, check=True)
                self.logger.info("Created BLAST database for prokka_variants")
                
    def run_prokka_variants_search(self):
        """Run BLAST search for prokka variants."""
        self.logger.info("Running prokka variants BLAST search...")
        self.prepare_blast_db()
        
        batch_genomes = self.get_batch_genomes(use_assemblies=False)
        
        for genome_dir in batch_genomes:
            if not genome_dir.is_dir():
                continue
                
            prefix = self.detect_prefix(genome_dir)
            if not prefix:
                self.logger.warning(f"Skipping {genome_dir}: no valid files found")
                continue
                
            genome_name = prefix
            ffn_files = list(genome_dir.glob("*.ffn"))
            ffn_file = ffn_files[0] if ffn_files else None
            
            if ffn_file and ffn_file.exists():
                self.run_blastn_variants(genome_name, ffn_file)
    
    def run_tblastn(self, genome_name: str, assembly_file: Path):
        """Run tblastn: protein query vs genome nucleotide."""
        output_file = self.blast_dirs['coding_protein'] / f"{genome_name}_genes_blast.txt"
        
        if output_file.exists() and output_file.stat().st_size > 0:
            self.logger.debug(f"Skipping tblastn for {genome_name}: output exists")
            return
            
        # Handle compressed assembly files
        if assembly_file.suffix == '.gz':
            # Use process substitution for compressed files
            import tempfile
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
                subprocess.run(['gunzip', '-c', str(assembly_file)], stdout=open(tmp.name, 'w'), check=True)
                temp_file = Path(tmp.name)
            
            cmd = [
                'tblastn',
                '-query', str(self.query_files['protein']),
                '-subject', str(temp_file),
                '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs',
                '-evalue', '1e-5',
                '-num_threads', str(self.num_threads),
                '-out', str(output_file)
            ]
            
            try:
                subprocess.run(cmd, check=True)
            finally:
                temp_file.unlink()
        else:
            cmd = [
                'tblastn',
                '-query', str(self.query_files['protein']),
                '-subject', str(assembly_file),
                '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs',
                '-evalue', '1e-5',
                '-num_threads', str(self.num_threads),
                '-out', str(output_file)
            ]
            subprocess.run(cmd, check=True)
        
        self.logger.debug(f"Completed tblastn for {genome_name}")
        
    def run_blastn_coding(self, genome_name: str, assembly_file: Path):
        """Run blastn: nucleotide query vs genome nucleotide."""
        output_file = self.blast_dirs['coding_nt'] / f"{genome_name}_genes_blast.txt"
        
        if output_file.exists() and output_file.stat().st_size > 0:
            self.logger.debug(f"Skipping blastn_coding for {genome_name}: output exists")
            return
            
        # Handle compressed assembly files
        if assembly_file.suffix == '.gz':
            import tempfile
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
                subprocess.run(['gunzip', '-c', str(assembly_file)], stdout=open(tmp.name, 'w'), check=True)
                temp_file = Path(tmp.name)
            
            cmd = [
                'blastn',
                '-query', str(self.query_files['nt']),
                '-subject', str(temp_file),
                '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs',
                '-evalue', '1e-5',
                '-num_threads', str(self.num_threads),
                '-out', str(output_file)
            ]
            
            try:
                subprocess.run(cmd, check=True)
            finally:
                temp_file.unlink()
        else:
            cmd = [
                'blastn',
                '-query', str(self.query_files['nt']),
                '-subject', str(assembly_file),
                '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs',
                '-evalue', '1e-5',
                '-num_threads', str(self.num_threads),
                '-out', str(output_file)
            ]
            
            self.logger.debug(f"Running blastn_coding for {genome_name}")
            subprocess.run(cmd, check=True)
        
    def run_blastn_variants(self, genome_name: str, ffn_file: Path):
        """Run blastn: Prokka genes vs reference operon database."""
        output_file = self.blast_dirs['prokka_variants'] / f"{genome_name}_blast.txt"
        
        if output_file.exists() and output_file.stat().st_size > 0:
            self.logger.debug(f"Skipping blastn_variants for {genome_name}: output exists")
            return
            
        ref_db = self.output_dir / 'reference_db' / 'operon_genes_nt'
        cmd = [
            'blastn',
            '-query', str(ffn_file),
            '-db', str(ref_db),
            '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq',
            '-perc_identity', '80',
            '-num_threads', str(self.num_threads),
            '-out', str(output_file)
        ]
        
        self.logger.debug(f"Running blastn_variants for {genome_name}")
        subprocess.run(cmd, check=True)
        
    def run_blastn_noncoding(self, genome_name: str, assembly_file: Path):
        """Run blastn: noncoding/promoter query vs genome."""
        output_file = self.blast_dirs['noncoding'] / f"{genome_name}_noncoding_blast.txt"
        
        if output_file.exists() and output_file.stat().st_size > 0:
            self.logger.debug(f"Skipping blastn_noncoding for {genome_name}: output exists")
            return
            
        # Handle compressed assembly files
        if assembly_file.suffix == '.gz':
            import tempfile
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
                subprocess.run(['gunzip', '-c', str(assembly_file)], stdout=open(tmp.name, 'w'), check=True)
                temp_file = Path(tmp.name)
            
            cmd = [
                'blastn',
                '-query', str(self.query_files['noncoding']),
                '-subject', str(temp_file),
                '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs',
                '-evalue', '1e-5',
                '-num_threads', str(self.num_threads),
                '-out', str(output_file)
            ]
            
            self.logger.debug(f"Running blastn_noncoding for {genome_name}")
            try:
                subprocess.run(cmd, check=True)
            finally:
                temp_file.unlink()
        else:
            cmd = [
                'blastn',
                '-query', str(self.query_files['noncoding']),
                '-subject', str(assembly_file),
                '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs',
                '-evalue', '1e-5',
                '-num_threads', str(self.num_threads),
                '-out', str(output_file)
            ]
            
            self.logger.debug(f"Running blastn_noncoding for {genome_name}")
            subprocess.run(cmd, check=True)
        
    def process_blast_results(self):
        """Process BLAST results to identify genomes with complete operons."""
        self.logger.info("Processing BLAST results...")
        
        all_results = []
        
        # Get all unique genome IDs from all BLAST modes
        genome_ids = set()
        
        # Check coding_protein directory
        blast_dir = self.blast_dirs['coding_protein']
        gene_files = list(blast_dir.glob('*_genes_blast.txt'))
        noncoding_files = list(blast_dir.glob('*_noncoding_blast.txt'))
        for f in gene_files:
            genome_id = f.stem.replace('_genes_blast', '')
            genome_ids.add(genome_id)
        for f in noncoding_files:
            genome_id = f.stem.replace('_noncoding_blast', '')
            genome_ids.add(genome_id)
            
        # Check coding_nt directory
        nt_dir = self.blast_dirs['coding_nt']
        if nt_dir.exists():
            for f in nt_dir.glob('*_genes_blast.txt'):
                genome_id = f.stem.replace('_genes_blast', '')
                genome_ids.add(genome_id)
                
        # Check prokka_variants directory
        variants_dir = self.blast_dirs['prokka_variants']
        if variants_dir.exists():
            for f in variants_dir.glob('*_blast.txt'):
                genome_id = f.stem.replace('_blast', '')
                genome_ids.add(genome_id)
            
        genome_ids = sorted(list(genome_ids))
        
            
        self.logger.info(f"Found {len(genome_ids)} genomes to process")
        
        # Track high-quality hits statistics
        hq_stats = {'total_hits': 0, 'hq_hits': 0}
        
        for genome_id in genome_ids:
            # Process tblastn results (coding_protein)
            gene_file = blast_dir / f'{genome_id}_genes_blast.txt'
            if gene_file.exists() and gene_file.stat().st_size > 0:
                df = pd.read_csv(gene_file, sep='\t', header=None,
                               names=['qseqid', 'sseqid', 'pident', 'length', 
                                     'mismatch', 'gapopen', 'qstart', 'qend',
                                     'sstart', 'send', 'evalue', 'bitscore', 'qcovs'])
                # Keep best hit per query gene
                best_hits = df.sort_values('bitscore', ascending=False).groupby('qseqid').first()
                best_hits['genome_id'] = genome_id
                best_hits['search_type'] = 'tblastn'
                # Add HQ flag
                best_hits['is_hq'] = (best_hits['pident'] >= 90) & (best_hits['qcovs'] >= 80)
                hq_stats['total_hits'] += len(best_hits)
                hq_stats['hq_hits'] += best_hits['is_hq'].sum()
                all_results.append(best_hits)
                
            # Process blastn results (coding_nt)
            nt_file = nt_dir / f'{genome_id}_genes_blast.txt'
            if nt_file.exists() and nt_file.stat().st_size > 0:
                df = pd.read_csv(nt_file, sep='\t', header=None,
                               names=['qseqid', 'sseqid', 'pident', 'length', 
                                     'mismatch', 'gapopen', 'qstart', 'qend',
                                     'sstart', 'send', 'evalue', 'bitscore', 'qcovs'])
                best_hits = df.sort_values('bitscore', ascending=False).groupby('qseqid').first()
                best_hits['genome_id'] = genome_id
                best_hits['search_type'] = 'blastn_nt'
                best_hits['is_hq'] = (best_hits['pident'] >= 90) & (best_hits['qcovs'] >= 80)
                hq_stats['total_hits'] += len(best_hits)
                hq_stats['hq_hits'] += best_hits['is_hq'].sum()
                all_results.append(best_hits)
                
            # Process non-coding BLAST results
            noncoding_file = blast_dir / f'{genome_id}_noncoding_blast.txt'
            if noncoding_file.exists() and noncoding_file.stat().st_size > 0:
                df = pd.read_csv(noncoding_file, sep='\t', header=None,
                               names=['qseqid', 'sseqid', 'pident', 'length', 
                                     'mismatch', 'gapopen', 'qstart', 'qend',
                                     'sstart', 'send', 'evalue', 'bitscore', 'qcovs'])
                best_hits = df.sort_values('bitscore', ascending=False).groupby('qseqid').first()
                best_hits['genome_id'] = genome_id
                best_hits['search_type'] = 'blastn_noncoding'
                best_hits['is_hq'] = (best_hits['pident'] >= 90) & (best_hits['qcovs'] >= 80)
                hq_stats['total_hits'] += len(best_hits)
                hq_stats['hq_hits'] += best_hits['is_hq'].sum()
                all_results.append(best_hits)
                
        if all_results:
            combined_df = pd.concat(all_results)
            combined_df.to_csv(self.output_dir / 'all_blast_hits.csv')
            self.logger.info(f"Saved combined results to {self.output_dir / 'all_blast_hits.csv'}")
            
            # Log HQ statistics
            if hq_stats['total_hits'] > 0:
                hq_pct = (hq_stats['hq_hits'] / hq_stats['total_hits']) * 100
                self.logger.info(f"High-quality hits (≥90% identity, ≥80% coverage): {hq_stats['hq_hits']}/{hq_stats['total_hits']} ({hq_pct:.1f}%)")
            
            # Generate summary
            self.generate_summary(combined_df)
            
        return combined_df if all_results else pd.DataFrame()
        
    def generate_summary(self, df: pd.DataFrame):
        """Generate summary statistics from BLAST results."""
        self.logger.info("Generating summary statistics...")
        
        # Count hits per genome and gene
        summary = []
        for genome_id in df['genome_id'].unique():
            genome_df = df[df['genome_id'] == genome_id]
            # Count hits from different search types
            gene_hits = genome_df[genome_df['search_type'].isin(['tblastn', 'blastn_nt'])]['qseqid'].nunique()
            noncoding_hits = genome_df[genome_df['search_type'] == 'blastn_noncoding']['qseqid'].nunique()
            # Count HQ hits
            hq_gene_hits = genome_df[(genome_df['search_type'].isin(['tblastn', 'blastn_nt'])) & 
                                    (genome_df['is_hq'] == True)]['qseqid'].nunique()
            
            summary.append({
                'genome_id': genome_id,
                'gene_hits': gene_hits,
                'gene_hits_hq': hq_gene_hits,
                'noncoding_hits': noncoding_hits,
                'total_hits': gene_hits + noncoding_hits
            })
            
        summary_df = pd.DataFrame(summary)
        summary_df = summary_df.sort_values('total_hits', ascending=False)
        summary_df.to_csv(self.output_dir / 'operon_summary.csv', index=False)
        
        # Print statistics
        self.logger.info("\n=== BLAST Summary Statistics ===")
        self.logger.info(f"Total genomes analyzed: {len(summary_df)}")
        self.logger.info(f"Genomes with complete operon (7 genes): {len(summary_df[summary_df['gene_hits'] == 7])}")
        self.logger.info(f"Genomes with complete HQ operon (7 genes, ≥90% identity, ≥80% coverage): {len(summary_df[summary_df['gene_hits_hq'] == 7])}")
        self.logger.info(f"Genomes with ≥6 genes: {len(summary_df[summary_df['gene_hits'] >= 6])}")
        self.logger.info(f"Genomes with ≥5 genes: {len(summary_df[summary_df['gene_hits'] >= 5])}")
        self.logger.info(f"Genomes with promoter hits: {len(summary_df[summary_df['noncoding_hits'] > 0])}")
        
    def create_blast_overview(self):
        """Create overview of BLAST results showing prevalence of each operon component."""
        self.logger.info("Creating BLAST overview...")
        
        all_hits = []
        genome_ids = set()
        
        # Get all unique genome IDs from all BLAST modes
        blast_dir = self.blast_dirs['coding_protein']
        gene_files = list(blast_dir.glob('*_genes_blast.txt'))
        noncoding_files = list(blast_dir.glob('*_noncoding_blast.txt'))
        for f in gene_files:
            genome_id = f.stem.replace('_genes_blast', '')
            genome_ids.add(genome_id)
            
        # Check coding_nt directory
        nt_dir = self.blast_dirs['coding_nt']
        if nt_dir.exists():
            for f in nt_dir.glob('*_genes_blast.txt'):
                genome_id = f.stem.replace('_genes_blast', '')
                genome_ids.add(genome_id)
                
        genome_ids = sorted(list(genome_ids))
        
            
        self.logger.info(f"Processing {len(genome_ids)} genomes...")
        
        for genome_id in genome_ids:
            # Process tblastn results
            gene_file = blast_dir / f'{genome_id}_genes_blast.txt'
            if gene_file.exists() and gene_file.stat().st_size > 0:
                df = pd.read_csv(gene_file, sep='\t', header=None,
                               names=['qseqid', 'sseqid', 'pident', 'length', 
                                     'mismatch', 'gapopen', 'qstart', 'qend',
                                     'sstart', 'send', 'evalue', 'bitscore', 'qcovs'])
                df['genome_id'] = genome_id
                df['search_type'] = 'tblastn'
                df['is_hq'] = (df['pident'] >= 90) & (df['qcovs'] >= 80)
                all_hits.append(df)
                
            # Process blastn_nt results
            nt_file = nt_dir / f'{genome_id}_genes_blast.txt'
            if nt_file.exists() and nt_file.stat().st_size > 0:
                df = pd.read_csv(nt_file, sep='\t', header=None,
                               names=['qseqid', 'sseqid', 'pident', 'length', 
                                     'mismatch', 'gapopen', 'qstart', 'qend',
                                     'sstart', 'send', 'evalue', 'bitscore', 'qcovs'])
                df['genome_id'] = genome_id
                df['search_type'] = 'blastn_nt'
                df['is_hq'] = (df['pident'] >= 90) & (df['qcovs'] >= 80)
                all_hits.append(df)
                
            # Process non-coding BLAST results
            noncoding_file = blast_dir / f'{genome_id}_noncoding_blast.txt'
            if noncoding_file.exists() and noncoding_file.stat().st_size > 0:
                df = pd.read_csv(noncoding_file, sep='\t', header=None,
                               names=['qseqid', 'sseqid', 'pident', 'length', 
                                     'mismatch', 'gapopen', 'qstart', 'qend',
                                     'sstart', 'send', 'evalue', 'bitscore', 'qcovs'])
                df['genome_id'] = genome_id
                df['search_type'] = 'blastn_noncoding'
                df['is_hq'] = (df['pident'] >= 90) & (df['qcovs'] >= 80)
                all_hits.append(df)
                
        if not all_hits:
            self.logger.warning("No BLAST hits found")
            return pd.DataFrame()
            
        # Combine all hits
        all_hits_df = pd.concat(all_hits, ignore_index=True)
        
        # Gene prevalence analysis (combine tblastn and blastn_nt results)
        coding_hits = all_hits_df[all_hits_df['search_type'].isin(['tblastn', 'blastn_nt'])]
        gene_prevalence = coding_hits.groupby('qseqid')['genome_id'].nunique()
        
        # Also calculate HQ prevalence
        hq_coding_hits = coding_hits[coding_hits['is_hq'] == True]
        gene_prevalence_hq = hq_coding_hits.groupby('qseqid')['genome_id'].nunique()
        
        # Align HQ counts with gene prevalence index
        hq_counts = gene_prevalence_hq.reindex(gene_prevalence.index).fillna(0).astype(int)
        
        gene_prevalence_df = pd.DataFrame({
            'gene': gene_prevalence.index,
            'genomes_with_hit': gene_prevalence.values,
            'genomes_with_hq_hit': hq_counts.values,
            'prevalence_percent': (gene_prevalence.values / len(genome_ids)) * 100,
            'hq_prevalence_percent': (hq_counts.values / len(genome_ids)) * 100
        })
        gene_prevalence_df = gene_prevalence_df.sort_values('prevalence_percent', ascending=False)
        
        # Save results
        gene_prevalence_df.to_csv(self.output_dir / 'gene_prevalence.csv', index=False)
        all_hits_df.to_csv(self.output_dir / 'all_blast_hits_detailed.csv', index=False)
        
        self.logger.info("\n=== Gene Prevalence ===")
        for _, row in gene_prevalence_df.iterrows():
            self.logger.info(f"{row['gene']}: {row['genomes_with_hit']}/{len(genome_ids)} ({row['prevalence_percent']:.1f}%) | HQ: {row['genomes_with_hq_hit']}/{len(genome_ids)} ({row['hq_prevalence_percent']:.1f}%)")
            
        return all_hits_df
        
    def generate_manuscript_stats(self):
        """Generate comprehensive statistics for manuscript."""
        self.logger.info("Generating manuscript statistics...")
        
        stats_file = self.output_dir / 'manuscript_stats.txt'
        
        with open(stats_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("BLAST SEARCH MANUSCRIPT STATISTICS\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("=" * 80 + "\n\n")
            
            # Check reference sequences
            f.write("REFERENCE SEQUENCES\n")
            f.write("-" * 40 + "\n")
            for name, path in self.query_files.items():
                if path.exists():
                    seq_count = sum(1 for line in open(path) if line.startswith('>'))
                    f.write(f"{name}: {seq_count} sequences\n")
                else:
                    f.write(f"{name}: File not found\n")
            f.write("\n")
            
            # BLAST output statistics
            f.write("BLAST OUTPUT FILES\n")
            f.write("-" * 40 + "\n")
            for mode, dir_path in self.blast_dirs.items():
                if dir_path.exists():
                    if mode == 'coding_protein':
                        file_count = len(list(dir_path.glob('*_genes_blast.txt')))
                    elif mode == 'noncoding':
                        file_count = len(list(dir_path.glob('*_noncoding_blast.txt')))
                    else:
                        file_count = len(list(dir_path.glob('*.txt')))
                    f.write(f"{mode}: {file_count} result files\n")
                else:
                    f.write(f"{mode}: Directory not found\n")
            f.write("\n")
            
            # Load and analyze results
            if (self.output_dir / 'operon_summary.csv').exists():
                summary_df = pd.read_csv(self.output_dir / 'operon_summary.csv')
                
                f.write("OPERON DETECTION SUMMARY\n")
                f.write("-" * 40 + "\n")
                f.write(f"Total genomes analyzed: {len(summary_df)}\n")
                f.write(f"Genomes with complete operon (7 genes): {len(summary_df[summary_df['gene_hits'] == 7])}\n")
                if 'gene_hits_hq' in summary_df.columns:
                    f.write(f"Genomes with complete HQ operon (7 genes, ≥90% identity, ≥80% coverage): {len(summary_df[summary_df['gene_hits_hq'] == 7])}\n")
                f.write(f"Genomes with ≥6 genes: {len(summary_df[summary_df['gene_hits'] >= 6])}\n")
                f.write(f"Genomes with ≥5 genes: {len(summary_df[summary_df['gene_hits'] >= 5])}\n")
                f.write(f"Genomes with any hits: {len(summary_df[summary_df['total_hits'] > 0])}\n")
                f.write("\n")
                
                # Distribution of gene counts
                f.write("GENE COUNT DISTRIBUTION\n")
                f.write("-" * 40 + "\n")
                for i in range(8):
                    count = len(summary_df[summary_df['gene_hits'] == i])
                    percent = (count / len(summary_df)) * 100
                    f.write(f"{i} genes: {count} genomes ({percent:.1f}%)\n")
                    
        self.logger.info(f"Manuscript statistics saved to {stats_file}")
        
def main():
    parser = argparse.ArgumentParser(description='Unified BLAST pipeline for operon analysis')
    parser.add_argument('--mode', choices=['search', 'process', 'overview', 'stats', 'all'],
                       default='all', help='Pipeline mode to run')
    parser.add_argument('--blast-modes', nargs='+', 
                       choices=['coding_protein', 'coding_nt', 'prokka_variants', 'noncoding'],
                       help='BLAST search modes to run')
    parser.add_argument('--batch-id', type=int, help='Batch ID for SLURM array job')
    parser.add_argument('--batch-size', type=int, default=100, help='Number of genomes per batch')
    parser.add_argument('--threads', type=int, default=60, help='Number of threads for BLAST')
    parser.add_argument('--prokka-dir', default='../01_prokka_annotation/output/prokka_results',
                       help='Directory containing Prokka results')
    parser.add_argument('--assembly-dir', default='../../Efs_assemblies',
                       help='Directory containing genome assemblies')
    parser.add_argument('--ref-dir', default='../02_reference_operon_extraction/output',
                       help='Directory containing reference sequences')
    parser.add_argument('--output-dir', default='output', help='Output directory')
    
    args = parser.parse_args()
    
    # Initialize pipeline
    pipeline = BlastPipeline(
        prokka_dir=args.prokka_dir,
        assembly_dir=args.assembly_dir,
        reference_dir=args.ref_dir,
        output_dir=args.output_dir,
        batch_id=args.batch_id,
        batch_size=args.batch_size,
        num_threads=args.threads
    )
    
    # Create output directories
    pipeline.create_directories()
    
    # Run requested mode(s)
    if args.mode in ['search', 'all']:
        pipeline.run_blast_search(modes=args.blast_modes)
        
    if args.mode in ['process', 'all']:
        pipeline.process_blast_results()
        
    if args.mode in ['overview', 'all']:
        pipeline.create_blast_overview()
        
    if args.mode in ['stats', 'all']:
        pipeline.generate_manuscript_stats()
        
    pipeline.logger.info("Pipeline completed successfully")
    
if __name__ == '__main__':
    main()