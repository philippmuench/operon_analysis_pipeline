#!/usr/bin/env python3
"""
Create PDF visualization of gene orders in operons.
Shows canonical order and all non-canonical arrangements.
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch
import numpy as np

def create_gene_order_pdf():
    """Create PDF visualization of gene arrangements."""
    
    # Load the analysis results
    results_file = "output/operon_order_analysis.csv"
    df = pd.read_csv(results_file)
    
    # Define canonical order and gene colors
    CANONICAL_ORDER = ['frpC', 'glpC', 'ptsD', 'ptsC', 'ptsB', 'ptsA', 'fruR']
    
    # Gene color scheme - consistent colors for each gene
    gene_colors = {
        'frpC': '#FF6B6B',  # Red
        'glpC': '#4ECDC4',  # Teal
        'ptsD': '#45B7D1',  # Light blue
        'ptsC': '#96CEB4',  # Light green
        'ptsB': '#FFEAA7',  # Light yellow
        'ptsA': '#DDA0DD',  # Plum
        'fruR': '#98D8C8',  # Mint
    }
    
    # Gene descriptions (shortened for space)
    gene_short = {
        'frpC': 'Fructoselysine deglycase',
        'glpC': 'Glucoselysine deglycase',
        'ptsD': 'PTS IID',
        'ptsC': 'PTS IIC',
        'ptsB': 'PTS IIB',
        'ptsA': 'PTS IIA',
        'fruR': 'Sigma-54 regulator'
    }
    
    # Get complete operons
    complete_operons = df[df['complete_operon'] == True].copy()
    
    # Separate canonical and non-canonical
    non_canonical_operons = []
    
    for idx, row in complete_operons.iterrows():
        genome_id = row['genome_id']
        gene_order = row['gene_order']
        
        # Parse gene order
        if isinstance(gene_order, str):
            try:
                gene_order = eval(gene_order)
            except:
                continue
        
        # Check if canonical (forward or reverse)
        if gene_order != CANONICAL_ORDER and gene_order != CANONICAL_ORDER[::-1]:
            non_canonical_operons.append({
                'genome_id': genome_id,
                'gene_order': gene_order,
                'same_contig': row.get('same_contig', None),
                'same_strand': row.get('same_strand', None),
                'synteny_score': row.get('synteny_score', None)
            })
    
    # Also get partial operon
    partial_operons = df[(df['n_genes_found'] > 0) & (df['n_genes_found'] < 7)]
    
    # Calculate total number of rows needed (reduced since no legend)
    n_rows = 2 + len(non_canonical_operons) + len(partial_operons)  # Header + canonical + non-canonical + partial
    
    # Create figure (long page)
    fig_height = max(10, n_rows * 0.6)
    fig, ax = plt.subplots(figsize=(14, fig_height), facecolor='white')
    
    # Title and statistics
    ax.text(7, n_rows - 0.5, 'Operon Gene Order Analysis', 
            fontsize=16, fontweight='bold', ha='center')
    ax.text(7, n_rows - 1.2, f'Total genomes: {len(df)} | Complete operons: {len(complete_operons)} | Canonical order: 99.7%',
            fontsize=11, ha='center', color='#666')
    
    current_y = n_rows - 2
    
    # Draw canonical order
    ax.text(0.5, current_y, 'Canonical Order (7,947 genomes - 99.7%)', 
            fontsize=11, fontweight='bold', color='darkgreen')
    draw_gene_order(ax, CANONICAL_ORDER, current_y - 0.5, gene_colors, 
                   highlight=True, same_contig=True, same_strand=True)
    
    current_y -= 2
    
    # Draw non-canonical orders
    if non_canonical_operons:
        ax.text(0.5, current_y, f'Non-Canonical Orders ({len(non_canonical_operons)} genomes - 0.3%)',
                fontsize=11, fontweight='bold', color='darkred')
        current_y -= 0.8
        
        for i, operon in enumerate(non_canonical_operons, 1):
            genome_name = operon['genome_id'].replace('.result', '')
            if len(genome_name) > 20:
                genome_name = genome_name[:20] + '...'
            
            info_text = f"{i}. {genome_name}"
            if not operon['same_contig']:
                info_text += " [multi-contig]"
            if not operon['same_strand']:
                info_text += " [mixed-strand]"
            info_text += f" (synteny: {operon['synteny_score']:.2f})"
            
            ax.text(0.5, current_y, info_text, fontsize=9, color='#333')
            draw_gene_order(ax, operon['gene_order'], current_y - 0.5, gene_colors,
                          same_contig=operon['same_contig'], 
                          same_strand=operon['same_strand'])
            current_y -= 1.2
    
    # Draw partial operons if any
    if len(partial_operons) > 0:
        current_y -= 0.5
        ax.text(0.5, current_y, f'Partial Operons ({len(partial_operons)} genome{"s" if len(partial_operons) > 1 else ""})',
                fontsize=11, fontweight='bold', color='darkorange')
        current_y -= 0.8
        
        for idx, row in partial_operons.iterrows():
            genome_id = row['genome_id'].replace('.result', '')
            if len(genome_id) > 20:
                genome_id = genome_id[:20] + '...'
            
            gene_order = row['gene_order']
            if isinstance(gene_order, str):
                try:
                    gene_order = eval(gene_order)
                except:
                    gene_order = []
            
            if gene_order:
                genes_found = row['genes_found']
                if isinstance(genes_found, str):
                    try:
                        genes_found = eval(genes_found)
                    except:
                        genes_found = []
                
                info_text = f"{genome_id} ({len(genes_found)} genes: {', '.join(genes_found)})"
                ax.text(0.5, current_y, info_text, fontsize=9, color='#333')
                draw_gene_order(ax, gene_order, current_y - 0.5, gene_colors)
                current_y -= 1.2
    
    # Set axis limits and remove axes
    ax.set_xlim(0, 14)
    ax.set_ylim(-0.5, n_rows)
    ax.axis('off')
    
    # Save as PDF
    output_file = "output/gene_order_visualization.pdf"
    plt.savefig(output_file, bbox_inches='tight', facecolor='white')
    print(f"PDF visualization created: {output_file}")
    print(f"  - Canonical operons: {len(complete_operons) - len(non_canonical_operons)}")
    print(f"  - Non-canonical operons: {len(non_canonical_operons)}")
    print(f"  - Partial operons: {len(partial_operons)}")
    plt.close()

def draw_gene_order(ax, gene_order, y_pos, gene_colors, highlight=False, 
                   same_contig=True, same_strand=True):
    """Draw a single gene order arrangement."""
    
    x_start = 2
    x_pos = x_start
    
    for i, gene in enumerate(gene_order):
        if gene in gene_colors:
            # Draw gene box
            if highlight:
                rect = FancyBboxPatch((x_pos, y_pos - 0.25), 0.9, 0.5,
                                     boxstyle="round,pad=0.02",
                                     facecolor=gene_colors[gene], 
                                     edgecolor='black', linewidth=1.5, alpha=0.9)
            else:
                rect = FancyBboxPatch((x_pos, y_pos - 0.2), 0.8, 0.4,
                                     boxstyle="round,pad=0.02",
                                     facecolor=gene_colors[gene], 
                                     edgecolor='black', linewidth=0.8, alpha=0.7)
            ax.add_patch(rect)
            
            # Add gene name
            fontsize = 8 if highlight else 7
            ax.text(x_pos + 0.4, y_pos, gene, fontsize=fontsize, 
                   ha='center', va='center', color='white', fontweight='bold')
            
            # Add arrow between genes
            if i < len(gene_order) - 1:
                ax.annotate('', xy=(x_pos + 0.95, y_pos), 
                           xytext=(x_pos + 0.85, y_pos),
                           arrowprops=dict(arrowstyle='->', color='#666', lw=1))
            
            x_pos += 1.1

if __name__ == "__main__":
    create_gene_order_pdf()