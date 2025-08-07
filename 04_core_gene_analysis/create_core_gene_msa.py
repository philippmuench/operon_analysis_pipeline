#!/usr/bin/env python3
"""
Create multiple sequence alignments for core genes using MAFFT.
"""

import os
import glob
import subprocess
from multiprocessing import Pool, cpu_count
from Bio import AlignIO, SeqIO
import argparse

def create_msa_for_gene(args):
    """Create MSA for a single gene using MAFFT."""
    gene_file, output_dir, timeout = args
    
    gene_name = os.path.basename(gene_file).replace('.fasta', '')
    output_file = os.path.join(output_dir, f"{gene_name}_aligned.fasta")
    
    # Check if input file has sequences
    try:
        sequences = list(SeqIO.parse(gene_file, "fasta"))
        if len(sequences) < 2:
            return gene_name, f"Skipped (only {len(sequences)} sequences)"
    except:
        return gene_name, "Error reading sequences"
    
    # Run MAFFT
    try:
        cmd = ["mafft", "--auto", "--quiet", gene_file]
        # Handle timeout=0 as no timeout (None)
        actual_timeout = None if timeout == 0 else timeout
        with open(output_file, 'w') as outf:
            result = subprocess.run(cmd, stdout=outf, stderr=subprocess.PIPE, 
                                    text=True, timeout=actual_timeout)
        
        if result.returncode == 0:
            # Verify alignment was created properly
            try:
                alignment = AlignIO.read(output_file, "fasta")
                return gene_name, f"Success ({len(alignment)} sequences, {alignment.get_alignment_length()} bp)"
            except:
                return gene_name, "Error: Invalid alignment produced"
        else:
            return gene_name, f"MAFFT failed: {result.stderr}"
            
    except subprocess.TimeoutExpired:
        return gene_name, "Timeout (>5 minutes)"
    except FileNotFoundError:
        return gene_name, "Error: MAFFT not found"
    except Exception as e:
        return gene_name, f"Error: {str(e)}"

def main():
    parser = argparse.ArgumentParser(description="Create MSAs for core genes")
    parser.add_argument("--sequences-dir", default="output/core_gene_sequences",
                        help="Directory containing core gene sequences")
    parser.add_argument("--output-dir", default="output/core_gene_alignments",
                        help="Output directory for alignments")
    parser.add_argument("--threads", type=int, default=min(cpu_count()-1, 10),
                        help="Number of threads to use")
    parser.add_argument("--timeout", type=int, default=None,
                        help="Timeout for MAFFT in seconds (default: no timeout)")
    
    args = parser.parse_args()
    
    # Check if MAFFT_TMPDIR is set (strict validation)
    if "MAFFT_TMPDIR" not in os.environ:
        print("❌ Error: MAFFT_TMPDIR environment variable is not set!")
        print("   MAFFT may use /tmp which can cause disk space issues on HPC systems.")
        print("   Please set it via: export MAFFT_TMPDIR=/vol/tmp")
        print("   Run this script from the main pipeline which sets this automatically.")
        return 1
    
    # Validate MAFFT_TMPDIR directory exists and is writable
    mafft_tmpdir = os.environ["MAFFT_TMPDIR"]
    if not os.path.exists(mafft_tmpdir):
        print(f"❌ Error: MAFFT_TMPDIR directory does not exist: {mafft_tmpdir}")
        print("   Please create the directory or set a different MAFFT_TMPDIR")
        return 1
    
    if not os.access(mafft_tmpdir, os.W_OK):
        print(f"❌ Error: MAFFT_TMPDIR is not writable: {mafft_tmpdir}")
        print("   Please check permissions on the temporary directory")
        return 1
    
    print(f"✅ MAFFT_TMPDIR validated: {mafft_tmpdir}")

    # Check if sequences directory exists
    if not os.path.exists(args.sequences_dir):
        print(f"Error: Sequences directory not found: {args.sequences_dir}")
        print("Run extract_core_sequences.py first!")
        return 1
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Find all gene sequence files
    gene_files = glob.glob(os.path.join(args.sequences_dir, "*.fasta"))
    
    if not gene_files:
        print(f"Error: No FASTA files found in {args.sequences_dir}")
        return 1
    
    print(f"{'='*60}")
    print(f"Creating MSAs for {len(gene_files)} core genes")
    print(f"{'='*60}")
    print(f"Input directory: {args.sequences_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Threads: {args.threads}")
    print(f"Timeout: {args.timeout or 'None'}")
    print(f"{'='*60}\n")
    
    # Prepare arguments for parallel processing
    task_args = [(gene_file, args.output_dir, args.timeout) for gene_file in gene_files]
    
    # Process in parallel
    print("Creating alignments...")
    with Pool(args.threads) as pool:
        results = pool.map(create_msa_for_gene, task_args)
    
    # Summarize results
    successful = 0
    failed = 0
    
    print(f"\n{'='*60}")
    print("MSA Creation Results:")
    print(f"{'='*60}")
    
    for gene_name, status in sorted(results):
        if status.startswith("Success"):
            successful += 1
            print(f"✅ {gene_name}: {status}")
        else:
            failed += 1
            print(f"❌ {gene_name}: {status}")
    
    print(f"\n{'='*60}")
    print(f"Summary:")
    print(f"  Successful: {successful}")
    print(f"  Failed: {failed}")
    print(f"  Total: {len(results)}")
    print(f"{'='*60}")
    
    # Create summary file
    summary_file = os.path.join(args.output_dir, "msa_summary.txt")
    with open(summary_file, 'w') as f:
        f.write("Core Gene MSA Creation Summary\n")
        f.write("==============================\n\n")
        f.write(f"Total genes processed: {len(results)}\n")
        f.write(f"Successful alignments: {successful}\n")
        f.write(f"Failed alignments: {failed}\n\n")
        f.write("Detailed results:\n")
        for gene_name, status in sorted(results):
            f.write(f"  {gene_name}: {status}\n")
    
    print(f"Summary saved to: {summary_file}")
    
    if successful > 0:
        print(f"\n✅ Successfully created {successful} MSAs!")
        return 0
    else:
        print(f"\n❌ No successful alignments created")
        return 1

if __name__ == "__main__":
    exit(main())