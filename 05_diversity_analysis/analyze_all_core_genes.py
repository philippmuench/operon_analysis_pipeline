#!/usr/bin/env python3
"""
Analyze diversity for ALL core genes (not just a subset).
Fixed version that doesn't require seqkit.
"""

import os
import pandas as pd
from Bio import SeqIO, AlignIO
import subprocess
from multiprocessing import Pool, cpu_count
import random
import numpy as np

def sample_sequences(input_file, output_file, n_sample=1000):
    """Sample sequences using Python instead of seqkit"""
    sequences = list(SeqIO.parse(input_file, 'fasta'))
    n_seqs = len(sequences)
    
    if n_seqs <= n_sample:
        # If we have fewer sequences than requested, just copy the file
        with open(output_file, 'w') as out:
            SeqIO.write(sequences, out, 'fasta')
    else:
        # Sample randomly
        sampled = random.sample(sequences, n_sample)
        with open(output_file, 'w') as out:
            SeqIO.write(sampled, out, 'fasta')

def calculate_diversity_single_gene(args):
    """Calculate diversity for a single gene"""
    gene_file, gene_name = args
    
    try:
        # Count sequences
        n_seqs = sum(1 for _ in SeqIO.parse(gene_file, 'fasta'))
        
        if n_seqs < 2:
            return None
            
        # For very large alignments, subsample
        if n_seqs > 1000:
            # Subsample to 1000 sequences
            subsample_file = f"{gene_file}.subsample"
            sample_sequences(gene_file, subsample_file, 1000)
            input_file = subsample_file
            subsampled = True
        else:
            input_file = gene_file
            subsampled = False
            
        # Run MAFFT with fast settings
        aligned_file = f"{gene_file}.aligned"
        cmd = f"mafft --retree 1 --maxiterate 0 {input_file} > {aligned_file} 2>/dev/null"
        subprocess.run(cmd, shell=True, check=True)
        
        # Calculate pairwise identity
        alignment = AlignIO.read(aligned_file, 'fasta')
        identities = []
        
        # Sample pairs for large alignments
        n_align_seqs = len(alignment)
        if n_align_seqs > 100:
            # Sample 5000 random pairs
            pairs = []
            for _ in range(5000):
                i = random.randint(0, n_align_seqs-1)
                j = random.randint(0, n_align_seqs-1)
                if i != j:
                    pairs.append((i, j))
        else:
            # All pairs
            pairs = [(i, j) for i in range(n_align_seqs) for j in range(i+1, n_align_seqs)]
        
        for i, j in pairs:
            seq1 = str(alignment[i].seq).upper()
            seq2 = str(alignment[j].seq).upper()
            
            # Skip gaps
            matches = sum(1 for a, b in zip(seq1, seq2) if a == b and a != '-')
            total = sum(1 for a, b in zip(seq1, seq2) if a != '-' and b != '-')
            
            if total > 0:
                identities.append(matches / total * 100)
        
        # Cleanup
        os.remove(aligned_file)
        if subsampled:
            os.remove(subsample_file)
        
        if identities:
            return {
                'gene': gene_name,
                'n_sequences': n_seqs,
                'avg_identity': np.mean(identities),
                'std_identity': np.std(identities),
                'min_identity': np.min(identities),
                'max_identity': np.max(identities),
                'n_comparisons': len(identities)
            }
    
    except Exception as e:
        print(f"Error processing {gene_name}: {str(e)}")
        return None

def main():
    print("=== Analyzing ALL Core Genes ===")
    print(f"Time: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Read core genes list
    with open('core_genes_95pct.txt', 'r') as f:
        core_genes = [line.strip() for line in f if line.strip()]
    
    print(f"Total core genes to analyze: {len(core_genes)}")
    
    # Get all sequence files
    seq_dir = 'core_gene_sequences_95pct'
    gene_files = []
    
    for gene in core_genes:
        seq_file = os.path.join(seq_dir, f"{gene.replace('/', '_').replace(' ', '_')}.fasta")
        if os.path.exists(seq_file):
            gene_files.append((seq_file, gene))
    
    print(f"Found sequence files for {len(gene_files)} genes")
    
    # Process in parallel
    n_cores = min(cpu_count() - 1, 80)
    print(f"Processing with {n_cores} cores...")
    
    results = []
    batch_size = 50
    n_batches = (len(gene_files) + batch_size - 1) // batch_size
    
    for batch_idx in range(n_batches):
        start_idx = batch_idx * batch_size
        end_idx = min((batch_idx + 1) * batch_size, len(gene_files))
        batch = gene_files[start_idx:end_idx]
        
        print(f"\nProcessing batch {batch_idx+1}/{n_batches} ({len(batch)} genes)...")
        
        with Pool(n_cores) as pool:
            batch_results = pool.map(calculate_diversity_single_gene, batch)
            results.extend([r for r in batch_results if r is not None])
        
        print(f"  Completed: {len(results)} genes so far")
    
    # Save results
    df = pd.DataFrame(results)
    df = df.sort_values('avg_identity', ascending=False)
    df.to_csv('all_core_genes_diversity.csv', index=False)
    
    print(f"\nAnalyzed {len(df)} genes successfully")
    print(f"Average identity: {df['avg_identity'].mean():.2f}% ± {df['avg_identity'].std():.2f}%")
    print(f"Range: {df['avg_identity'].min():.2f}% - {df['avg_identity'].max():.2f}%")
    
    # Compare with operon genes
    operon_genes = {
        "fructoselysine-6-phosphate_deglycase": 99.60,
        "glucoselysine-6-phosphate_deglycase": 98.96,
        "fructoselysine_glucoselysine_PTS_system_EIID_component": 99.38,
        "fructoselysine_glucoselysine_PTS_system_EIIC_component": 99.26,
        "fructoselysine_glucoselysine_PTS_system_EIIB_component_[EC:2.7.1.-]": 99.42,
        "fructoselysine_glucoselysine_PTS_system_EIIA_component_[EC:2.7.1.-]": 99.51,
        "sigma-54_dependent_transcriptional_regulator,_gfr_operon_transcriptional_activator": 98.69
    }
    
    print("\nOperon gene rankings among ALL core genes:")
    for gene, identity in operon_genes.items():
        # Find percentile
        percentile = (df['avg_identity'] < identity).sum() / len(df) * 100
        print(f"  {gene}: {identity:.2f}% ({percentile:.1f} percentile)")
    
    # Create visualization
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot core gene distribution
    ax.hist(df['avg_identity'], bins=50, alpha=0.7, color='blue', label='Core genes')
    
    # Add operon genes as vertical lines
    colors = plt.cm.Set1(range(len(operon_genes)))
    for i, (gene, identity) in enumerate(operon_genes.items()):
        short_name = gene.split('_')[0] if '_' in gene else gene[:20]
        ax.axvline(identity, color=colors[i], linestyle='--', linewidth=2, 
                   label=f'{short_name} ({identity:.1f}%)')
    
    ax.set_xlabel('Average Pairwise Identity (%)', fontsize=12)
    ax.set_ylabel('Number of Genes', fontsize=12)
    ax.set_title(f'Conservation of ALL {len(df)} Core Genes vs Operon Genes', fontsize=14)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig('all_core_vs_operon_diversity.png', dpi=300, bbox_inches='tight')
    print("\nVisualization saved to all_core_vs_operon_diversity.png")

if __name__ == '__main__':
    main()