#!/usr/bin/env python3
"""
Unified Prokka Pipeline Management Tool
Combines all functionality for running, checking, validating, and re-running Prokka annotations.
"""

import os
import sys
import glob
import gzip
import subprocess
import argparse
import logging
from datetime import datetime
from pathlib import Path
import json
import numpy as np

class ProkkaPipeline:
    def __init__(self, test_mode=False):
        self.test_mode = test_mode
        self.script_dir = Path(__file__).parent.absolute()
        self.genome_dir = self.script_dir.parent.parent / "Efs_assemblies"
        
        if test_mode:
            self.output_dir = self.script_dir / "output" / "prokka_results"
            self.genome_list_file = self.script_dir / "genome_list_test.txt"
        else:
            self.output_dir = self.script_dir.parent / "prokka_output"
            self.genome_list_file = self.script_dir / "genome_list.txt"
        
        self.expected_files = ["gff", "faa", "ffn", "fna", "gbk", "log", "txt", "tsv"]
        
        # Setup logging
        log_file = self.script_dir / f"prokka_pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)
    
    def create_genome_list(self):
        """Create list of genome files."""
        genome_files = sorted(glob.glob(str(self.genome_dir / "*.fasta.gz")))
        
        if self.test_mode:
            # In test mode, only use first 50 genomes
            genome_files = genome_files[:50]
            output_file = self.genome_list_file
        else:
            output_file = self.script_dir / "genome_list.txt"
        
        with open(output_file, 'w') as f:
            for genome in genome_files:
                f.write(f"{genome}\n")
        
        self.logger.info(f"Created genome list with {len(genome_files)} genomes: {output_file}")
        return len(genome_files)
    
    def run_prokka_batch(self, start_idx, end_idx):
        """Run Prokka on a batch of genomes."""
        # Read genome list
        with open(self.genome_list_file, 'r') as f:
            all_genomes = [line.strip() for line in f if line.strip()]
        
        # Get batch
        batch_genomes = all_genomes[start_idx:end_idx]
        
        self.logger.info(f"Processing batch: genomes {start_idx+1} to {min(end_idx, len(all_genomes))}")
        
        success_count = 0
        failed_genomes = []
        
        for genome_path in batch_genomes:
            basename = Path(genome_path).stem.replace('.fasta', '')
            outdir = self.output_dir / basename
            
            # Skip if already processed
            if (outdir / f"{basename}.gff").exists():
                self.logger.info(f"Skipping {basename} - already processed")
                success_count += 1
                continue
            
            # Extract genome
            temp_fasta = self.script_dir / f"temp_{os.getpid()}_{basename}.fasta"
            try:
                with gzip.open(genome_path, 'rt') as f_in:
                    with open(temp_fasta, 'w') as f_out:
                        f_out.write(f_in.read())
            except Exception as e:
                self.logger.error(f"Failed to extract {basename}: {e}")
                failed_genomes.append(genome_path)
                continue
            
            # Run Prokka
            cmd = [
                "prokka",
                "--outdir", str(outdir),
                "--prefix", basename,
                "--kingdom", "Bacteria",
                "--genus", "Enterococcus",
                "--species", "faecalis",
                "--cpus", str(os.environ.get('SLURM_CPUS_PER_TASK', 4)),
                "--centre", "X",
                "--compliant",
                "--force",
                "--quiet",
                str(temp_fasta)
            ]
            
            try:
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
                if result.returncode == 0:
                    self.logger.info(f"Successfully processed {basename}")
                    success_count += 1
                else:
                    self.logger.error(f"Prokka failed for {basename}: {result.stderr}")
                    failed_genomes.append(genome_path)
            except subprocess.TimeoutExpired:
                self.logger.error(f"Prokka timeout for {basename}")
                failed_genomes.append(genome_path)
            except Exception as e:
                self.logger.error(f"Error running Prokka for {basename}: {e}")
                failed_genomes.append(genome_path)
            finally:
                # Clean up temp file
                if temp_fasta.exists():
                    temp_fasta.unlink()
        
        return success_count, failed_genomes
    
    def check_progress(self):
        """Check Prokka annotation progress."""
        # In test mode, count from genome list file if it exists
        if self.test_mode and self.genome_list_file.exists():
            with open(self.genome_list_file, 'r') as f:
                total_genomes = sum(1 for line in f if line.strip())
        else:
            total_genomes = len(list(self.genome_dir.glob("*.fasta.gz")))
        
        if self.output_dir.exists():
            completed = len(list(self.output_dir.glob("*/*.gff")))
        else:
            completed = 0
        
        remaining = total_genomes - completed
        progress = (completed / total_genomes * 100) if total_genomes > 0 else 0
        
        results = {
            "total_genomes": total_genomes,
            "completed": completed,
            "remaining": remaining,
            "progress_percentage": progress
        }
        
        self.logger.info("Prokka Annotation Progress")
        self.logger.info("=" * 50)
        self.logger.info(f"Total genomes: {total_genomes}")
        self.logger.info(f"Completed: {completed}")
        self.logger.info(f"Remaining: {remaining}")
        self.logger.info(f"Progress: {progress:.1f}%")
        
        return results
    
    def validate_outputs(self):
        """Validate Prokka annotation outputs."""
        # Get genome list
        if self.genome_list_file.exists():
            with open(self.genome_list_file, 'r') as f:
                genome_list = [line.strip() for line in f if line.strip()]
        else:
            genome_list = sorted(glob.glob(str(self.genome_dir / "*.fasta.gz")))
        
        total_genomes = len(genome_list)
        complete_genomes = 0
        incomplete_genomes = []
        missing_outputs = 0
        empty_files = 0
        
        self.logger.info("Prokka Output Validation Report")
        self.logger.info("=" * 50)
        self.logger.info(f"Total genomes to validate: {total_genomes}")
        
        for genome_path in genome_list:
            basename = Path(genome_path).stem.replace('.fasta', '')
            genome_dir = self.output_dir / basename
            
            genome_complete = True
            genome_issues = []
            
            # Check if output directory exists
            if not genome_dir.exists():
                self.logger.warning(f"Missing output directory for {basename}")
                incomplete_genomes.append(genome_path)
                missing_outputs += len(self.expected_files)
                continue
            
            # Check each expected file
            for ext in self.expected_files:
                file_path = genome_dir / f"{basename}.{ext}"
                
                if not file_path.exists():
                    genome_complete = False
                    genome_issues.append(f"Missing: {basename}.{ext}")
                    missing_outputs += 1
                elif file_path.stat().st_size == 0:
                    genome_complete = False
                    genome_issues.append(f"Empty: {basename}.{ext}")
                    empty_files += 1
            
            if genome_complete:
                complete_genomes += 1
            else:
                incomplete_genomes.append(genome_path)
                if genome_issues:
                    self.logger.warning(f"INCOMPLETE: {basename}")
                    for issue in genome_issues:
                        self.logger.warning(f"  - {issue}")
        
        # Summary statistics
        self.logger.info("")
        self.logger.info("Summary Statistics")
        self.logger.info("=" * 50)
        self.logger.info(f"Total genomes: {total_genomes}")
        self.logger.info(f"Complete genomes: {complete_genomes} ({complete_genomes/total_genomes*100:.1f}%)")
        self.logger.info(f"Incomplete genomes: {len(incomplete_genomes)}")
        self.logger.info(f"Missing files: {missing_outputs}")
        self.logger.info(f"Empty files: {empty_files}")
        
        # Save list of incomplete genomes
        if incomplete_genomes:
            incomplete_list_file = self.script_dir / f"incomplete_genomes_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
            with open(incomplete_list_file, 'w') as f:
                for genome in incomplete_genomes:
                    f.write(f"{genome}\n")
            self.logger.info(f"Incomplete genomes list saved to: {incomplete_list_file}")
        
        return {
            "total_genomes": total_genomes,
            "complete_genomes": complete_genomes,
            "incomplete_genomes": len(incomplete_genomes),
            "incomplete_list": incomplete_genomes,
            "missing_files": missing_outputs,
            "empty_files": empty_files
        }
    
    def rerun_failed(self, incomplete_list_file=None):
        """Rerun Prokka for failed/incomplete genomes."""
        # Find incomplete list file
        if not incomplete_list_file:
            # Find most recent incomplete genomes list
            incomplete_files = sorted(self.script_dir.glob("incomplete_genomes_*.txt"))
            if not incomplete_files:
                self.logger.error("No incomplete genomes list found. Run validation first.")
                return None
            incomplete_list_file = incomplete_files[-1]
        
        incomplete_list_file = Path(incomplete_list_file)
        if not incomplete_list_file.exists():
            self.logger.error(f"Incomplete list file not found: {incomplete_list_file}")
            return None
        
        # Read list of incomplete genomes
        with open(incomplete_list_file, 'r') as f:
            incomplete_genomes = [line.strip() for line in f if line.strip()]
        
        self.logger.info(f"Reprocessing {len(incomplete_genomes)} failed/incomplete genomes")
        
        success_count = 0
        still_failed = []
        
        for genome_path in incomplete_genomes:
            if not Path(genome_path).exists():
                self.logger.warning(f"Genome file not found: {genome_path}")
                still_failed.append(genome_path)
                continue
            
            basename = Path(genome_path).stem.replace('.fasta', '')
            outdir = self.output_dir / basename
            
            # Remove any existing incomplete output
            if outdir.exists():
                self.logger.info(f"Removing incomplete output for {basename}")
                import shutil
                shutil.rmtree(outdir)
            
            # Extract and run Prokka
            temp_fasta = self.script_dir / f"temp_rerun_{os.getpid()}_{basename}.fasta"
            try:
                with gzip.open(genome_path, 'rt') as f_in:
                    with open(temp_fasta, 'w') as f_out:
                        f_out.write(f_in.read())
                
                cmd = [
                    "prokka",
                    "--outdir", str(outdir),
                    "--prefix", basename,
                    "--kingdom", "Bacteria",
                    "--genus", "Enterococcus",
                    "--species", "faecalis",
                    "--cpus", str(os.environ.get('SLURM_CPUS_PER_TASK', 4)),
                    "--force",
                    "--quiet",
                    str(temp_fasta)
                ]
                
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
                if result.returncode == 0:
                    self.logger.info(f"Successfully reprocessed {basename}")
                    success_count += 1
                else:
                    self.logger.error(f"Prokka failed again for {basename}")
                    still_failed.append(genome_path)
            except Exception as e:
                self.logger.error(f"Error reprocessing {basename}: {e}")
                still_failed.append(genome_path)
            finally:
                if temp_fasta.exists():
                    temp_fasta.unlink()
        
        self.logger.info("")
        self.logger.info("Reprocessing Summary")
        self.logger.info("=" * 50)
        self.logger.info(f"Total genomes processed: {len(incomplete_genomes)}")
        self.logger.info(f"Successful: {success_count}")
        self.logger.info(f"Still failed: {len(still_failed)}")
        
        return {
            "total_processed": len(incomplete_genomes),
            "successful": success_count,
            "still_failed": still_failed
        }
    
    def generate_manuscript_stats(self):
        """Generate manuscript statistics."""
        stats = {}
        
        # Count input genomes
        input_count = len(list(self.genome_dir.glob("*.fasta.gz")))
        stats['input_count'] = input_count
        
        # Get Prokka version
        try:
            result = subprocess.run(['prokka', '--version'], 
                                  capture_output=True, text=True, timeout=10)
            version = result.stderr.strip() if result.stderr else result.stdout.strip()
            stats['prokka_version'] = version
        except:
            stats['prokka_version'] = "Unknown"
        
        # Analyze outputs
        if not self.output_dir.exists():
            stats['error'] = "Output directory not found"
            return stats
        
        genome_dirs = [d for d in self.output_dir.iterdir() if d.is_dir()]
        stats['total_genomes'] = len(genome_dirs)
        
        complete_genomes = 0
        total_genes = 0
        gene_counts = []
        
        for genome_dir in genome_dirs:
            # Check if genome is complete
            complete = True
            for ext in self.expected_files:
                file_path = genome_dir / f"{genome_dir.name}.{ext}"
                if not file_path.exists() or file_path.stat().st_size == 0:
                    complete = False
                    break
            
            if complete:
                complete_genomes += 1
                
                # Count genes from GFF file
                gff_file = genome_dir / f"{genome_dir.name}.gff"
                if gff_file.exists():
                    try:
                        with open(gff_file, 'r') as f:
                            gene_count = sum(1 for line in f if not line.startswith('#') and 'CDS' in line)
                        gene_counts.append(gene_count)
                        total_genes += gene_count
                    except:
                        pass
        
        stats['complete_genomes'] = complete_genomes
        stats['total_genes'] = total_genes
        
        if gene_counts:
            stats['avg_genes_per_genome'] = sum(gene_counts) / len(gene_counts)
            stats['median_genes_per_genome'] = np.median(gene_counts)
            stats['std_genes_per_genome'] = np.std(gene_counts)
            stats['min_genes_per_genome'] = min(gene_counts)
            stats['max_genes_per_genome'] = max(gene_counts)
        
        # Print summary
        self.logger.info("")
        self.logger.info("MANUSCRIPT STATISTICS")
        self.logger.info("=" * 50)
        self.logger.info(f"Input assemblies: {input_count:,} E. faecalis genomes")
        self.logger.info(f"Prokka version: {stats['prokka_version']}")
        self.logger.info(f"Successfully annotated: {complete_genomes:,} genomes")
        if gene_counts:
            self.logger.info(f"Average genes per genome: {stats['avg_genes_per_genome']:.0f}")
            self.logger.info(f"Total genes predicted: {total_genes:,}")
        
        # Save to file
        stats_file = self.script_dir / f"manuscript_stats_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=2)
        self.logger.info(f"Statistics saved to: {stats_file}")
        
        return stats

def main():
    parser = argparse.ArgumentParser(description='Unified Prokka Pipeline Management Tool')
    parser.add_argument('command', choices=['run', 'check', 'validate', 'rerun', 'stats', 'prepare'],
                       help='Command to execute')
    parser.add_argument('--test', action='store_true', help='Run in test mode (first 50 genomes)')
    parser.add_argument('--batch-start', type=int, help='Start index for batch processing (0-based)')
    parser.add_argument('--batch-end', type=int, help='End index for batch processing')
    parser.add_argument('--incomplete-list', help='Path to incomplete genomes list file')
    parser.add_argument('--array-task-id', type=int, help='SLURM array task ID')
    parser.add_argument('--array-task-max', type=int, help='SLURM array task max')
    
    args = parser.parse_args()
    
    # Initialize pipeline
    pipeline = ProkkaPipeline(test_mode=args.test)
    
    if args.command == 'prepare':
        # Create genome list
        count = pipeline.create_genome_list()
        print(f"Created genome list with {count} genomes")
    
    elif args.command == 'run':
        # Handle SLURM array job or manual batch
        if args.array_task_id:
            # Calculate batch for SLURM array job
            batch_size = 100
            start_idx = (args.array_task_id - 1) * batch_size
            end_idx = args.array_task_id * batch_size
        elif args.batch_start is not None and args.batch_end is not None:
            start_idx = args.batch_start
            end_idx = args.batch_end
        else:
            print("Error: Must specify either --array-task-id or both --batch-start and --batch-end")
            sys.exit(1)
        
        success, failed = pipeline.run_prokka_batch(start_idx, end_idx)
        print(f"Batch completed: {success} successful, {len(failed)} failed")
        
        # Generate stats if this is the last array job
        if args.array_task_id and args.array_task_max:
            if args.array_task_id == args.array_task_max or args.test:
                pipeline.generate_manuscript_stats()
    
    elif args.command == 'check':
        results = pipeline.check_progress()
        print(json.dumps(results, indent=2))
    
    elif args.command == 'validate':
        results = pipeline.validate_outputs()
        print(f"Validation complete: {results['complete_genomes']}/{results['total_genomes']} genomes complete")
    
    elif args.command == 'rerun':
        results = pipeline.rerun_failed(args.incomplete_list)
        if results:
            print(f"Rerun complete: {results['successful']}/{results['total_processed']} successful")
    
    elif args.command == 'stats':
        stats = pipeline.generate_manuscript_stats()
        print(json.dumps(stats, indent=2))

if __name__ == "__main__":
    main()