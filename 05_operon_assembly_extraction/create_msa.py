#!/usr/bin/env python3
"""
Enhanced MSA creation for both coding and non-coding operon sequences.
Creates multiple sequence alignments using MAFFT.
"""

import os
import subprocess
import os
import pandas as pd
import argparse
from Bio import SeqIO, AlignIO

def run_mafft(input_file, output_file, threads=1):
    """Run MAFFT alignment."""
    cmd = [
        'mafft',
        '--auto',
        '--thread', str(threads),
        input_file
    ]
    
    try:
        with open(output_file, 'w') as f:
            env = os.environ.copy()
            # Default MAFFT tmp to /vol/tmp if not already set
            env.setdefault('MAFFT_TMPDIR', '/vol/tmp')
            env.setdefault('TMPDIR', env['MAFFT_TMPDIR'])
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True, env=env)
        
        if result.returncode != 0:
            print(f"MAFFT failed for {input_file}: {result.stderr}")
            return False
        return True
    except Exception as e:
        print(f"Error running MAFFT for {input_file}: {e}")
        return False

def create_alignments(sequences_dir, output_dir, threads=1, sequence_type="coding"):
    """Create multiple sequence alignments for all sequence files."""
    
    # Create output directory
    if sequence_type == "coding":
        align_dir = os.path.join(output_dir, "dna_alignments")
    else:
        align_dir = os.path.join(output_dir, "noncoding_alignments")
    
    os.makedirs(align_dir, exist_ok=True)
    
    # Get sequence files
    if os.path.exists(sequences_dir):
        seq_files = [f for f in os.listdir(sequences_dir) if f.endswith('.fasta')]
    else:
        print(f"Warning: Sequences directory {sequences_dir} not found")
        return []
    
    alignment_stats = []
    successful_alignments = []
    
    for seq_file in seq_files:
        gene_name = seq_file.replace('.fasta', '')
        input_path = os.path.join(sequences_dir, seq_file)
        output_path = os.path.join(align_dir, f"{gene_name}_aligned.fasta")
        
        print(f"Creating alignment for {gene_name}...")
        
        # Check if we have enough sequences
        seq_count = 0
        seq_lengths = []
        
        with open(input_path, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                seq_count += 1
                seq_lengths.append(len(record.seq))
        
        if seq_count < 2:
            print(f"Skipping {gene_name}: insufficient sequences ({seq_count})")
            continue
        
        # Run MAFFT
        if run_mafft(input_path, output_path, threads):
            # Calculate alignment statistics
            alignment = AlignIO.read(output_path, "fasta")
            alignment_length = alignment.get_alignment_length()
            n_sequences = len(alignment)
            
            # Calculate some basic stats
            avg_seq_length = sum(seq_lengths) / len(seq_lengths)
            min_seq_length = min(seq_lengths)
            max_seq_length = max(seq_lengths)
            
            stats = {
                'gene': gene_name,
                'type': sequence_type,
                'n_sequences': n_sequences,
                'alignment_length': alignment_length,
                'avg_original_length': round(avg_seq_length, 2),
                'min_original_length': min_seq_length,
                'max_original_length': max_seq_length,
                'alignment_file': output_path
            }
            
            alignment_stats.append(stats)
            successful_alignments.append(output_path)
            
            print(f"  Aligned {n_sequences} sequences, length: {alignment_length}")
        else:
            print(f"  Failed to align {gene_name}")
    
    return alignment_stats, successful_alignments

def main():
    parser = argparse.ArgumentParser(description="Create enhanced MSAs for coding and non-coding sequences")
    parser.add_argument("--coding-sequences", default="output/operon_sequences", 
                        help="Directory containing coding gene sequences")
    parser.add_argument("--noncoding-sequences", default="output/noncoding_sequences", 
                        help="Directory containing non-coding sequences")
    parser.add_argument("--output-dir", default="output/msa", 
                        help="Output directory for alignments")
    parser.add_argument("--threads", type=int, default=1, 
                        help="Number of threads for MAFFT")
    parser.add_argument("--coding-only", action="store_true", 
                        help="Only process coding sequences")
    parser.add_argument("--noncoding-only", action="store_true", 
                        help="Only process non-coding sequences")
    
    args = parser.parse_args()
    
    all_stats = []
    
    # Process coding sequences
    if not args.noncoding_only:
        print("Creating alignments for coding sequences...")
        coding_stats, coding_alignments = create_alignments(
            args.coding_sequences, 
            args.output_dir, 
            args.threads, 
            "coding"
        )
        all_stats.extend(coding_stats)
        print(f"Created {len(coding_alignments)} coding sequence alignments")
    
    # Process non-coding sequences  
    if not args.coding_only:
        print("\nCreating alignments for non-coding sequences...")
        noncoding_stats, noncoding_alignments = create_alignments(
            args.noncoding_sequences, 
            args.output_dir, 
            args.threads, 
            "noncoding"
        )
        all_stats.extend(noncoding_stats)
        print(f"Created {len(noncoding_alignments)} non-coding sequence alignments")
    
    # Save alignment summary
    if all_stats:
        summary_file = os.path.join(args.output_dir, "enhanced_alignment_summary.csv")
        summary_df = pd.DataFrame(all_stats)
        summary_df.to_csv(summary_file, index=False)
        print(f"\nAlignment summary saved to {summary_file}")
        
        # Print summary statistics
        print("\nAlignment Summary:")
        print(f"Total alignments created: {len(all_stats)}")
        
        for seq_type in ['coding', 'noncoding']:
            type_stats = summary_df[summary_df['type'] == seq_type]
            if not type_stats.empty:
                print(f"\n{seq_type.capitalize()} sequences:")
                print(f"  Number of alignments: {len(type_stats)}")
                print(f"  Average sequences per alignment: {type_stats['n_sequences'].mean():.1f}")
                print(f"  Average alignment length: {type_stats['alignment_length'].mean():.1f}")
    else:
        print("No alignments were created")

if __name__ == "__main__":
    main()