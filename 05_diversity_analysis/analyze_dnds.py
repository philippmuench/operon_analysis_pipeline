#!/usr/bin/env python3
"""
Calculate dN/dS ratios for genes to assess selection pressure.
"""

import os
import sys
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count

def get_codon_changes(codon1, codon2, genetic_code=11):
    """Classify codon changes as synonymous or non-synonymous."""
    try:
        aa1 = str(Seq(codon1).translate(table=genetic_code))
        aa2 = str(Seq(codon2).translate(table=genetic_code))
        
        if aa1 == aa2:
            return 'synonymous'
        elif aa1 == '*' or aa2 == '*':
            return 'stop'
        else:
            return 'non-synonymous'
    except:
        return 'error'

def calculate_dnds_for_alignment(alignment_file):
    """Calculate dN/dS for a single alignment."""
    try:
        # Read alignment
        alignment = AlignIO.read(alignment_file, 'fasta')
        gene_name = os.path.basename(alignment_file).replace('.msa', '').replace('.fasta', '')
        
        # Convert to codon alignment
        codon_positions = []
        for i in range(0, len(alignment[0].seq) - 2, 3):
            codon_positions.append(i)
        
        # Count substitutions
        synonymous = 0
        non_synonymous = 0
        pairs_analyzed = 0
        
        # Sample sequence pairs
        n_seqs = len(alignment)
        if n_seqs < 2:
            return None
            
        # For large alignments, sample pairs
        if n_seqs > 100:
            n_pairs = min(1000, n_seqs * (n_seqs - 1) // 2)
            pairs = []
            for _ in range(n_pairs):
                i = np.random.randint(0, n_seqs)
                j = np.random.randint(0, n_seqs)
                if i != j:
                    pairs.append((i, j))
        else:
            pairs = [(i, j) for i in range(n_seqs) for j in range(i+1, n_seqs)]
        
        for i, j in pairs:
            seq1 = str(alignment[i].seq).upper()
            seq2 = str(alignment[j].seq).upper()
            
            # Compare codons
            for pos in codon_positions:
                codon1 = seq1[pos:pos+3]
                codon2 = seq2[pos:pos+3]
                
                # Skip if gaps or incomplete codons
                if '-' in codon1 or '-' in codon2 or len(codon1) != 3 or len(codon2) != 3:
                    continue
                    
                if codon1 != codon2:
                    change_type = get_codon_changes(codon1, codon2)
                    if change_type == 'synonymous':
                        synonymous += 1
                    elif change_type == 'non-synonymous':
                        non_synonymous += 1
            
            pairs_analyzed += 1
        
        # Calculate dN/dS
        if synonymous > 0:
            dn_ds = non_synonymous / synonymous
        elif non_synonymous > 0:
            dn_ds = float('inf')
        else:
            dn_ds = 0
        
        return {
            'gene': gene_name,
            'n_sequences': n_seqs,
            'synonymous': synonymous,
            'non_synonymous': non_synonymous,
            'pairs_analyzed': pairs_analyzed,
            'dn_ds_ratio': dn_ds
        }
        
    except Exception as e:
        print(f"Error processing {alignment_file}: {str(e)}")
        return None

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Calculate dN/dS ratios')
    parser.add_argument('--input_dir', required=True, help='Directory with alignments')
    parser.add_argument('--output', required=True, help='Output CSV file')
    parser.add_argument('--cores', type=int, default=8, help='Number of cores')
    parser.add_argument('--file_pattern', default='*.msa', help='File pattern')
    
    args = parser.parse_args()
    
    # Get alignment files
    import glob
    alignment_files = glob.glob(os.path.join(args.input_dir, args.file_pattern))
    
    if not alignment_files:
        alignment_files = glob.glob(os.path.join(args.input_dir, '*.fasta'))
    
    print(f"Found {len(alignment_files)} alignment files")
    
    # Process in parallel
    with Pool(args.cores) as pool:
        results = pool.map(calculate_dnds_for_alignment, alignment_files)
    
    # Filter out None results
    results = [r for r in results if r is not None]
    
    # Create dataframe
    df = pd.DataFrame(results)
    df = df.sort_values('dn_ds_ratio')
    
    # Save results
    df.to_csv(args.output, index=False)
    
    # Print summary
    print(f"\nAnalyzed {len(df)} genes")
    print(f"Mean dN/dS: {df['dn_ds_ratio'][df['dn_ds_ratio'] != float('inf')].mean():.3f}")
    print(f"Genes under purifying selection (dN/dS < 1): {(df['dn_ds_ratio'] < 1).sum()}")
    print(f"Genes under positive selection (dN/dS > 1): {(df['dn_ds_ratio'] > 1).sum()}")
    
    # Selection categories
    strong_purifying = (df['dn_ds_ratio'] < 0.5).sum()
    weak_purifying = ((df['dn_ds_ratio'] >= 0.5) & (df['dn_ds_ratio'] < 1)).sum()
    neutral = ((df['dn_ds_ratio'] >= 0.9) & (df['dn_ds_ratio'] <= 1.1)).sum()
    positive = (df['dn_ds_ratio'] > 1.5).sum()
    
    print(f"\nSelection categories:")
    print(f"  Strong purifying (dN/dS < 0.5): {strong_purifying}")
    print(f"  Weak purifying (0.5 ≤ dN/dS < 1): {weak_purifying}")
    print(f"  Neutral (0.9 ≤ dN/dS ≤ 1.1): {neutral}")
    print(f"  Positive selection (dN/dS > 1.5): {positive}")

if __name__ == '__main__':
    main()