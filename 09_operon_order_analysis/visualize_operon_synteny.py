#!/usr/bin/env python3
"""
Create visualizations for operon order and synteny analysis.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
from matplotlib.patches import FancyBboxPatch, Rectangle
import json

# Set style for publication-quality figures
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

# Gene colors for consistency
GENE_COLORS = {
    'frpC': '#FF6B6B',  # Red
    'glpC': '#4ECDC4',  # Teal
    'ptsD': '#45B7D1',  # Blue
    'ptsC': '#96CEB4',  # Green
    'ptsB': '#FECA57',  # Yellow
    'ptsA': '#DDA0DD',  # Plum
    'fruR': '#98D8C8',  # Mint
}

CANONICAL_ORDER = ['frpC', 'glpC', 'ptsD', 'ptsC', 'ptsB', 'ptsA', 'fruR']

def create_synteny_diagram(df, output_dir, n_examples=10):
    """Create synteny diagram showing different operon organizations."""
    
    # Select representative examples of different operon types
    examples = []
    
    # Try to get examples of each type
    operon_types = ['Complete and conserved', 'Complete with rearrangement', 
                    'Complete but fragmented', 'Near-complete', 'Partial']
    
    for otype in operon_types:
        type_df = df[df['operon_type'] == otype]
        if len(type_df) > 0:
            # Get up to 2 examples of each type
            sample_size = min(2, len(type_df))
            examples.extend(type_df.sample(sample_size).to_dict('records'))
    
    # If we don't have enough examples, add more
    if len(examples) < n_examples:
        remaining = n_examples - len(examples)
        other_df = df[~df['genome_id'].isin([e['genome_id'] for e in examples])]
        if len(other_df) > 0:
            examples.extend(other_df.sample(min(remaining, len(other_df))).to_dict('records'))
    
    if not examples:
        print("No examples found for synteny diagram")
        return
    
    # Create figure
    fig, ax = plt.subplots(figsize=(14, len(examples) * 0.8), facecolor='white')
    
    # Draw canonical reference at top
    y_pos = len(examples) + 0.5
    draw_operon(ax, CANONICAL_ORDER, y_pos, label="Reference (Canonical)", is_reference=True)
    
    # Draw examples
    for i, example in enumerate(examples):
        y_pos = len(examples) - i - 0.5
        
        # Parse gene details if available
        gene_order = example['gene_order']
        if isinstance(gene_order, str):
            gene_order = eval(gene_order)  # Convert string representation to list
        
        label = f"{example['genome_id'][:15]}... ({example['operon_type']})"
        draw_operon(ax, gene_order, y_pos, label=label, 
                   same_strand=example.get('same_strand', True),
                   same_contig=example.get('same_contig', True))
    
    # Set labels and limits
    ax.set_xlim(-1, 9)
    ax.set_ylim(-0.5, len(examples) + 1)
    ax.set_xlabel('Gene Position', fontsize=11)
    ax.set_title('Operon Synteny Variations', fontsize=13, fontweight='bold')
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    # Add legend
    legend_elements = []
    for gene in CANONICAL_ORDER:
        legend_elements.append(patches.Patch(color=GENE_COLORS[gene], label=gene))
    ax.legend(handles=legend_elements, loc='upper right', ncol=7, fontsize=9,
             title="Genes", frameon=True, fancybox=True, shadow=True)
    
    plt.tight_layout()
    
    # Save figure
    output_file = os.path.join(output_dir, 'operon_synteny_examples.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(os.path.join(output_dir, 'operon_synteny_examples.pdf'), 
                bbox_inches='tight', facecolor='white')
    print(f"Saved synteny diagram to {output_file}")
    plt.close()

def draw_operon(ax, gene_order, y_pos, label="", same_strand=True, same_contig=True, is_reference=False):
    """Draw a single operon representation."""
    
    # Add label
    ax.text(-0.5, y_pos, label, fontsize=9, va='center', ha='right',
           fontweight='bold' if is_reference else 'normal')
    
    # Draw genes
    x_pos = 0
    for i, gene in enumerate(gene_order):
        if gene in GENE_COLORS:
            color = GENE_COLORS[gene]
            
            # Create gene box
            if is_reference:
                box = FancyBboxPatch((x_pos, y_pos - 0.3), 0.9, 0.6,
                                     boxstyle="round,pad=0.05",
                                     facecolor=color, edgecolor='black',
                                     linewidth=2, alpha=0.9)
            else:
                box = FancyBboxPatch((x_pos, y_pos - 0.25), 0.9, 0.5,
                                     boxstyle="round,pad=0.05",
                                     facecolor=color, edgecolor='black',
                                     linewidth=1, alpha=0.7)
            ax.add_patch(box)
            
            # Add gene name
            ax.text(x_pos + 0.45, y_pos, gene, fontsize=7, va='center', ha='center',
                   fontweight='bold', color='white')
            
            # Add arrow for strand direction
            if not same_strand and i % 2 == 1:  # Simulate inversions
                ax.annotate('', xy=(x_pos + 0.1, y_pos - 0.4), 
                           xytext=(x_pos + 0.8, y_pos - 0.4),
                           arrowprops=dict(arrowstyle='<-', color='red', lw=1))
            
            x_pos += 1.1
            
            # Add gap if not on same contig
            if not same_contig and i == 3:  # Simulate contig break
                ax.text(x_pos - 0.05, y_pos, '//', fontsize=12, va='center', 
                       ha='center', color='red', fontweight='bold')
                x_pos += 0.3

def create_operon_type_plot(df, output_dir):
    """Create bar plot showing operon organization types distribution."""
    
    fig, ax = plt.subplots(figsize=(10, 6), facecolor='white')
    
    # Get operon type counts
    type_counts = df['operon_type'].value_counts()
    
    # Create bar plot
    bars = ax.bar(range(len(type_counts)), type_counts.values,
                   color=sns.color_palette('viridis', len(type_counts)))
    
    # Customize plot
    ax.set_xticks(range(len(type_counts)))
    ax.set_xticklabels(type_counts.index, rotation=45, ha='right')
    ax.set_ylabel('Number of Genomes', fontsize=12)
    ax.set_xlabel('Operon Organization Type', fontsize=12)
    ax.set_title('Distribution of Operon Organization Types', fontsize=14, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)
    
    # Add value labels on bars
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height)}\n({100*height/len(df):.1f}%)', 
                ha='center', va='bottom', fontsize=10)
    
    # Add total count
    ax.text(0.02, 0.98, f'Total genomes: {len(df):,}', 
            transform=ax.transAxes, fontsize=11, va='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    
    # Save figure
    output_file = os.path.join(output_dir, 'operon_organization_types.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(os.path.join(output_dir, 'operon_organization_types.pdf'),
                bbox_inches='tight', facecolor='white')
    print(f"Saved operon types plot to {output_file}")
    plt.close()

def main():
    """Main function to create visualizations."""
    
    output_dir = 'output'
    
    # Load analysis results
    results_file = os.path.join(output_dir, 'operon_order_analysis.csv')
    
    if not os.path.exists(results_file):
        print(f"Error: Analysis results not found at {results_file}")
        print("Run analyze_operon_order.py first!")
        return
    
    print("Loading analysis results...")
    df = pd.read_csv(results_file)
    print(f"Loaded data for {len(df)} genomes")
    
    # Create visualizations
    print("\nCreating operon type distribution plot...")
    create_operon_type_plot(df, output_dir)
    
    print("\nVisualization complete!")

if __name__ == "__main__":
    main()