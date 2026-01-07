#!/usr/bin/env python3
"""
Regenerate detailed conservation metrics and plots from filtered MSAs.
Run from 05_operon_assembly_extraction directory.
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import SeqIO
from collections import Counter
from datetime import datetime


def calculate_shannon_entropy(alignment_file):
    """Calculate Shannon entropy (conservation) per position."""
    sequences = list(SeqIO.parse(alignment_file, 'fasta'))
    if not sequences:
        return [], {}

    aln_length = len(sequences[0].seq)
    n_seqs = len(sequences)

    conservation_scores = []
    gap_counts = []

    for pos in range(aln_length):
        column = [str(seq.seq[pos]).upper() for seq in sequences]
        counts = Counter(column)
        total = sum(counts.values())
        gaps = counts.get('-', 0) + counts.get('.', 0)
        gap_counts.append(gaps)

        # Shannon entropy
        entropy = 0
        for count in counts.values():
            if count > 0:
                p = count / total
                entropy -= p * np.log2(p)

        # Normalize and convert to conservation (1 - normalized_entropy)
        max_entropy = np.log2(min(21, total))  # 21 for amino acids + gap
        norm_entropy = entropy / max_entropy if max_entropy > 0 else 0
        conservation = 1 - norm_entropy
        conservation_scores.append(conservation)

    metadata = {
        'num_sequences': n_seqs,
        'alignment_length': aln_length,
        'mean_conservation': np.mean(conservation_scores),
        'mean_gaps': np.mean(gap_counts),
    }

    return conservation_scores, metadata


def calculate_detailed_metrics(alignment_file):
    """Calculate detailed conservation metrics."""
    sequences = list(SeqIO.parse(alignment_file, 'fasta'))
    if not sequences:
        return None

    aln_length = len(sequences[0].seq)
    n_seqs = len(sequences)

    conservation_scores = []
    gap_counts = []

    for pos in range(aln_length):
        column = [str(seq.seq[pos]).upper() for seq in sequences]
        counts = Counter(column)
        total = sum(counts.values())
        gaps = counts.get('-', 0) + counts.get('.', 0)
        gap_counts.append(gaps)

        entropy = 0
        for count in counts.values():
            if count > 0:
                p = count / total
                entropy -= p * np.log2(p)

        max_entropy = np.log2(min(21, total))
        norm_entropy = entropy / max_entropy if max_entropy > 0 else 0
        conservation = 1 - norm_entropy
        conservation_scores.append(conservation)

    # Categorize positions
    highly_conserved = sum(1 for c in conservation_scores if c >= 0.95)
    moderately_conserved = sum(1 for c in conservation_scores if 0.7 <= c < 0.95)
    variable = sum(1 for c in conservation_scores if c < 0.7)

    return {
        'n_sequences': n_seqs,
        'alignment_length': aln_length,
        'mean_conservation': np.mean(conservation_scores),
        'median_conservation': np.median(conservation_scores),
        'std_conservation': np.std(conservation_scores),
        'highly_conserved_pct': 100 * highly_conserved / aln_length,
        'moderately_conserved_pct': 100 * moderately_conserved / aln_length,
        'variable_pct': 100 * variable / aln_length,
        'mean_gaps_per_position': np.mean(gap_counts),
        'max_gaps_per_position': max(gap_counts),
        'positions_with_gaps_pct': 100 * sum(1 for g in gap_counts if g > 0) / aln_length,
    }


def create_conservation_plot(alignment_file, output_path, gene_name):
    """Create Shannon entropy conservation plot."""
    conservation_scores, metadata = calculate_shannon_entropy(alignment_file)
    if not conservation_scores:
        return

    positions = list(range(1, len(conservation_scores) + 1))

    fig, ax = plt.subplots(figsize=(15, 6))
    ax.plot(positions, conservation_scores, color='blue', linewidth=1.5, alpha=0.8)
    ax.fill_between(positions, conservation_scores, alpha=0.3)

    ax.set_xlabel('Position', fontsize=12)
    ax.set_ylabel('Conservation Score', fontsize=12)
    ax.set_title(f'{gene_name} Conservation (Shannon Entropy)\n'
                 f'n={metadata["num_sequences"]:,} sequences, '
                 f'mean={metadata["mean_conservation"]:.3f}', fontsize=14)
    ax.set_ylim(0, 1.05)
    ax.set_xlim(1, len(conservation_scores))
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


def main():
    print("=" * 60)
    print("Regenerating Conservation Metrics and Plots")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    msa_dir = 'output/msa/dna_alignments'
    plots_dir = 'output/plots'
    os.makedirs(plots_dir, exist_ok=True)

    # Find all alignment files
    msa_files = [f for f in os.listdir(msa_dir) if f.endswith('.fasta')]
    print(f"\nFound {len(msa_files)} MSA files")

    # Calculate detailed metrics
    results = []
    for msa_file in sorted(msa_files):
        filepath = os.path.join(msa_dir, msa_file)
        gene_name = msa_file.replace('_dna.fasta', '').replace('.fasta', '').replace('_aligned', '')

        print(f"  Processing {gene_name}...")
        metrics = calculate_detailed_metrics(filepath)
        if metrics:
            metrics['gene'] = gene_name
            results.append(metrics)

            # Create plot
            plot_path = os.path.join(plots_dir, f"{gene_name}_shannon_entropy_conservation.png")
            create_conservation_plot(filepath, plot_path, gene_name)

    # Save detailed metrics
    df = pd.DataFrame(results)
    cols = ['gene', 'n_sequences', 'alignment_length', 'mean_conservation',
            'median_conservation', 'std_conservation', 'highly_conserved_pct',
            'moderately_conserved_pct', 'variable_pct', 'mean_gaps_per_position',
            'max_gaps_per_position', 'positions_with_gaps_pct']
    df = df[cols]
    df.to_csv('output/operon_conservation_metrics_detailed.csv', index=False)

    print(f"\n✅ Saved output/operon_conservation_metrics_detailed.csv")
    print(f"✅ Saved {len(results)} conservation plots to {plots_dir}/")

    # Summary
    print("\nSummary:")
    for _, row in df.iterrows():
        print(f"  {row['gene']}: {row['n_sequences']} seqs, conservation={row['mean_conservation']:.4f}")

    print("\n" + "=" * 60)
    print("DONE!")
    print("=" * 60)


if __name__ == '__main__':
    main()
