#!/usr/bin/env python3
"""Calculate diversity metrics for core genes and compare with operon"""

import os
import glob
import pandas as pd
import numpy as np
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from collections import defaultdict
import subprocess
from multiprocessing import Pool, cpu_count

def calculate_pairwise_identity(alignment):
    """Calculate average pairwise identity for an alignment"""
    n_seqs = len(alignment)
    if n_seqs < 2:
        return 100.0
    
    identities = []
    for i in range(n_seqs):
        for j in range(i+1, n_seqs):
            seq1 = str(alignment[i].seq)
            seq2 = str(alignment[j].seq)
            
            # Count identical positions
            identical = sum(1 for a, b in zip(seq1, seq2) 
                          if a == b and a != '-' and b != '-')
            
            # Count comparable positions
            comparable = sum(1 for a, b in zip(seq1, seq2) 
                           if a != '-' and b != '-')
            
            if comparable > 0:
                identity = (identical / comparable) * 100
                identities.append(identity)
    
    return np.mean(identities) if identities else 100.0

def align_and_calculate_diversity(fasta_file):
    """Align sequences and calculate diversity metrics"""
    
    gene_name = os.path.basename(fasta_file).replace('.fasta', '')
    
    # Check if file exists and has sequences
    if not os.path.exists(fasta_file):
        return None
    
    sequences = list(SeqIO.parse(fasta_file, 'fasta'))
    if len(sequences) < 2:
        return None
    
    # Create alignment using MAFFT
    aligned_file = fasta_file.replace('.fasta', '_aligned.fasta')
    
    try:
        # Run MAFFT
        cmd = f"mafft --quiet --auto {fasta_file} > {aligned_file}"
        subprocess.run(cmd, shell=True, check=True, capture_output=True)
        
        # Read alignment
        alignment = AlignIO.read(aligned_file, 'fasta')
        
        # Calculate metrics
        avg_identity = calculate_pairwise_identity(alignment)
        
        # Calculate other diversity metrics
        seq_lengths = [len(str(seq.seq).replace('-', '')) for seq in alignment]
        avg_length = np.mean(seq_lengths)
        
        # Clean up
        os.remove(aligned_file)
        
        return {
            'gene': gene_name,
            'n_sequences': len(sequences),
            'avg_identity': avg_identity,
            'avg_length': avg_length,
            'std_length': np.std(seq_lengths)
        }
        
    except Exception as e:
        print(f"Error processing {gene_name}: {e}")
        return None

def process_batch(fasta_files):
    """Process a batch of fasta files"""
    results = []
    for fasta_file in fasta_files:
        result = align_and_calculate_diversity(fasta_file)
        if result:
            results.append(result)
    return results

def main():
    """Calculate diversity for core genes and compare with operon"""
    
    print("Calculating diversity metrics for core genes...")
    
    # Get all core gene sequence files
    fasta_files = glob.glob('core_gene_sequences/*.fasta')
    print(f"Found {len(fasta_files)} gene sequence files")
    
    # Sample if too many
    if len(fasta_files) > 200:
        import random
        random.seed(42)
        fasta_files = random.sample(fasta_files, 200)
        print(f"Sampling 200 genes for analysis")
    
    # Process in parallel
    n_cores = min(cpu_count() - 1, 8)
    batch_size = len(fasta_files) // n_cores + 1
    
    batches = [fasta_files[i:i+batch_size] for i in range(0, len(fasta_files), batch_size)]
    
    print(f"Processing with {n_cores} cores...")
    with Pool(n_cores) as pool:
        batch_results = pool.map(process_batch, batches)
    
    # Combine results
    all_results = []
    for batch in batch_results:
        all_results.extend(batch)
    
    # Create DataFrame
    diversity_df = pd.DataFrame(all_results)
    diversity_df = diversity_df.sort_values('avg_identity')
    
    # Save results
    diversity_df.to_csv('core_gene_diversity.csv', index=False)
    
    print(f"\nProcessed {len(diversity_df)} genes")
    print(f"Average identity across core genes: {diversity_df['avg_identity'].mean():.2f}%")
    print(f"Identity range: {diversity_df['avg_identity'].min():.2f}% - {diversity_df['avg_identity'].max():.2f}%")
    
    # Load operon diversity data
    print("\nComparing with operon genes...")
    operon_summary = pd.read_csv('operon_blast_summary.csv')
    
    # Calculate operon identities
    operon_genes = [
        'fructoselysine-6-phosphate_deglycase',
        'glucoselysine-6-phosphate_deglycase',
        'fructoselysine_glucoselysine_PTS_system_EIID_component',
        'fructoselysine_glucoselysine_PTS_system_EIIC_component',
        'fructoselysine_glucoselysine_PTS_system_EIIB_component_[EC:2.7.1.-]',
        'fructoselysine_glucoselysine_PTS_system_EIIA_component_[EC:2.7.1.-]',
        'sigma-54_dependent_transcriptional_regulator,_gfr_operon_transcriptional_activator'
    ]
    
    operon_identities = []
    for gene in operon_genes:
        col = f"{gene}_pident"
        if col in operon_summary.columns:
            # Get non-zero values
            values = operon_summary[col][operon_summary[col] > 0]
            if len(values) > 0:
                operon_identities.append({
                    'gene': gene.replace('_', ' '),
                    'avg_identity': values.mean(),
                    'n_sequences': len(values)
                })
    
    operon_df = pd.DataFrame(operon_identities)
    
    # Compare distributions
    print("\nOperon gene identities:")
    for _, row in operon_df.iterrows():
        print(f"  {row['gene']}: {row['avg_identity']:.2f}%")
    
    print(f"\nOperon average: {operon_df['avg_identity'].mean():.2f}%")
    print(f"Core genes average: {diversity_df['avg_identity'].mean():.2f}%")
    
    # Calculate percentile of operon genes in core gene distribution
    core_identities = diversity_df['avg_identity'].values
    
    print("\nOperon genes percentile in core gene distribution:")
    for _, operon_gene in operon_df.iterrows():
        percentile = np.sum(core_identities <= operon_gene['avg_identity']) / len(core_identities) * 100
        print(f"  {operon_gene['gene']}: {percentile:.1f} percentile")
    
    # Save comparison
    comparison_results = {
        'operon_mean_identity': operon_df['avg_identity'].mean(),
        'core_genes_mean_identity': diversity_df['avg_identity'].mean(),
        'core_genes_std_identity': diversity_df['avg_identity'].std(),
        'n_core_genes_analyzed': len(diversity_df),
        'operon_genes': operon_df.to_dict('records')
    }
    
    import json
    with open('diversity_comparison.json', 'w') as f:
        json.dump(comparison_results, f, indent=2)
    
    print("\nResults saved:")
    print("  - core_gene_diversity.csv")
    print("  - diversity_comparison.json")

if __name__ == "__main__":
    main()