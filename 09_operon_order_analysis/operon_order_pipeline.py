#!/usr/bin/env python3
"""
Unified pipeline for operon order and synteny analysis.
Consolidates all analysis steps into a single script.
Includes metadata stratification for analyzing patterns by source niche, continent, etc.
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, Rectangle
import seaborn as sns
import glob
from collections import defaultdict, Counter
import json
import argparse
from pathlib import Path

# Set style for publication-quality figures
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

# Define canonical operon gene order
CANONICAL_ORDER = ['frpC', 'glpC', 'ptsD', 'ptsC', 'ptsB', 'ptsA', 'fruR']

GENE_FUNCTIONS = {
    'frpC': 'Fructoselysine-6-phosphate deglycase',
    'glpC': 'Glucoselysine-6-phosphate deglycase',
    'ptsD': 'PTS system EIID component',
    'ptsC': 'PTS system EIIC component', 
    'ptsB': 'PTS system EIIB component',
    'ptsA': 'PTS system EIIA component',
    'fruR': 'Sigma-54 dependent transcriptional regulator'
}

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


class OperonOrderAnalyzer:
    """Main analyzer class for operon order and synteny analysis."""
    
    def __init__(self, blast_dir, output_dir, min_identity=80):
        self.blast_dir = blast_dir
        self.output_dir = output_dir
        self.min_identity = min_identity
        self.results_df = None
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
    
    def parse_blast_file(self, blast_file):
        """Parse a BLAST output file to extract hit information."""
        hits = []
        
        if not os.path.exists(blast_file):
            return hits
        
        try:
            # Read BLAST tabular output with 13 columns
            df = pd.read_csv(blast_file, sep='\t', header=None,
                            names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                                   'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                                   'evalue', 'bitscore', 'qcovs'])
            
            for _, row in df.iterrows():
                hit = {
                    'query': row['qseqid'],
                    'subject': row['sseqid'],
                    'identity': row['pident'],
                    'start': min(row['sstart'], row['send']),
                    'end': max(row['sstart'], row['send']),
                    'strand': '+' if row['sstart'] < row['send'] else '-',
                    'evalue': row['evalue'],
                    'bitscore': row['bitscore']
                }
                hits.append(hit)
        
        except Exception as e:
            print(f"Error parsing {blast_file}: {e}")
        
        return hits
    
    def extract_gene_from_query(self, query_id):
        """Extract gene name from query ID."""
        query_lower = query_id.lower()
        
        for gene in CANONICAL_ORDER:
            if gene.lower() in query_lower:
                return gene
        
        # Alternative names
        if 'fructoselysine' in query_lower and 'deglycase' in query_lower:
            return 'frpC'
        elif 'glucoselysine' in query_lower and 'deglycase' in query_lower:
            return 'glpC'
        elif 'pts' in query_lower:
            if 'eiid' in query_lower or 'iid' in query_lower:
                return 'ptsD'
            elif 'eiic' in query_lower or 'iic' in query_lower:
                return 'ptsC'
            elif 'eiib' in query_lower or 'iib' in query_lower:
                return 'ptsB'
            elif 'eiia' in query_lower or 'iia' in query_lower:
                return 'ptsA'
        elif 'sigma' in query_lower or 'regulator' in query_lower:
            return 'fruR'
        
        return None
    
    def calculate_synteny_score(self, observed_order):
        """Calculate synteny conservation score."""
        if not observed_order:
            return 0
        
        # Get indices of observed genes in canonical order
        canonical_indices = []
        for gene in observed_order:
            if gene in CANONICAL_ORDER:
                canonical_indices.append(CANONICAL_ORDER.index(gene))
        
        if len(canonical_indices) < 2:
            return 1.0 if len(canonical_indices) == 1 else 0
        
        # Check if order is preserved (forward or reverse)
        forward_score = 0
        reverse_score = 0
        
        for i in range(len(canonical_indices) - 1):
            # Forward direction
            if canonical_indices[i+1] > canonical_indices[i]:
                forward_score += 1
            # Reverse direction  
            if canonical_indices[i+1] < canonical_indices[i]:
                reverse_score += 1
        
        max_score = max(forward_score, reverse_score)
        return max_score / (len(canonical_indices) - 1)
    
    def analyze_genome_operon(self, genome_id, blast_files):
        """Analyze operon structure for a single genome."""
        
        # Collect all hits for operon genes
        gene_hits = defaultdict(list)
        
        for blast_file in blast_files:
            hits = self.parse_blast_file(blast_file)
            
            for hit in hits:
                gene = self.extract_gene_from_query(hit['query'])
                if gene and hit['identity'] >= self.min_identity and hit['bitscore'] >= 50:
                    gene_hits[gene].append(hit)
        
        if not gene_hits:
            return {
                'genome_id': genome_id,
                'genes_found': [],
                'n_genes_found': 0,
                'complete_operon': False,
                'gene_order': [],
                'canonical_order': False,
                'same_contig': False,
                'same_strand': False,
                'n_contigs': 0,
                'synteny_score': 0,
                'has_inversion': False,
                'has_rearrangement': False,
                'mean_gene_distance': None,
                'max_gene_distance': None,
                'operon_type': 'No genes',
                'gene_details': None
            }
        
        # Get best hit for each gene (highest bitscore)
        best_hits = {}
        for gene, hits in gene_hits.items():
            best_hit = max(hits, key=lambda x: x['bitscore'])
            best_hits[gene] = best_hit
        
        # Analyze operon structure
        genes_found = list(best_hits.keys())
        n_genes = len(genes_found)
        
        # Check completeness
        complete_operon = n_genes == len(CANONICAL_ORDER)
        
        # Analyze gene order
        gene_order = []
        gene_positions = []
        strands = []
        contigs = set()
        
        if n_genes > 0:
            # Sort genes by their genomic position
            gene_info = []
            for gene, hit in best_hits.items():
                gene_info.append({
                    'gene': gene,
                    'contig': hit['subject'],
                    'start': hit['start'],
                    'end': hit['end'],
                    'strand': hit['strand'],
                    'identity': hit['identity']
                })
                contigs.add(hit['subject'])
            
            # Sort by contig and position
            gene_info.sort(key=lambda x: (x['contig'], x['start']))
            
            gene_order = [g['gene'] for g in gene_info]
            gene_positions = [(g['start'], g['end']) for g in gene_info]
            strands = [g['strand'] for g in gene_info]
        
        # Check if genes are on same contig
        same_contig = len(contigs) == 1
        
        # Check if all genes are on same strand
        same_strand = len(set(strands)) == 1 if strands else False
        
        # Calculate synteny score
        synteny_score = self.calculate_synteny_score(gene_order)
        
        # Check for inversions or rearrangements
        has_inversion = False
        has_rearrangement = False
        
        if same_contig and n_genes > 1:
            # Check if order matches canonical
            canonical_indices = [CANONICAL_ORDER.index(g) for g in gene_order if g in CANONICAL_ORDER]
            if canonical_indices:
                diffs = np.diff(canonical_indices)
                if not all(d > 0 for d in diffs) and not all(d < 0 for d in diffs):
                    has_rearrangement = True
            
            # Check for strand switches
            if len(set(strands)) > 1:
                has_inversion = True
        
        # Calculate distances between adjacent genes
        gene_distances = []
        if same_contig and len(gene_positions) > 1:
            for i in range(len(gene_positions) - 1):
                dist = gene_positions[i+1][0] - gene_positions[i][1]
                gene_distances.append(dist)
        
        # Determine operon type
        if complete_operon and same_contig and same_strand and synteny_score > 0.9:
            operon_type = 'Complete and conserved'
        elif complete_operon and same_contig:
            operon_type = 'Complete with rearrangement'
        elif complete_operon:
            operon_type = 'Complete but fragmented'
        elif n_genes >= 5:
            operon_type = 'Near-complete'
        elif n_genes >= 3:
            operon_type = 'Partial'
        else:
            operon_type = 'Minimal'
        
        return {
            'genome_id': genome_id,
            'genes_found': genes_found,
            'n_genes_found': n_genes,
            'complete_operon': complete_operon,
            'gene_order': gene_order,
            'canonical_order': gene_order == CANONICAL_ORDER or gene_order == CANONICAL_ORDER[::-1],
            'same_contig': same_contig,
            'same_strand': same_strand,
            'n_contigs': len(contigs),
            'synteny_score': synteny_score,
            'has_inversion': has_inversion,
            'has_rearrangement': has_rearrangement,
            'mean_gene_distance': np.mean(gene_distances) if gene_distances else None,
            'max_gene_distance': max(gene_distances) if gene_distances else None,
            'operon_type': operon_type,
            'gene_details': json.dumps(gene_info) if 'gene_info' in locals() else None
        }
    
    def run_analysis(self):
        """Run the main operon order analysis."""
        
        print("="*60)
        print("Operon Order and Synteny Analysis")
        print("="*60)
        print(f"BLAST results directory: {self.blast_dir}")
        print(f"Canonical operon order: {' -> '.join(CANONICAL_ORDER)}")
        print()
        
        # Find BLAST files
        blast_files = glob.glob(os.path.join(self.blast_dir, "ENT_*_genes_blast.txt"))
        
        if not blast_files:
            print(f"Error: No BLAST files found in {self.blast_dir}")
            return False
        
        # Extract unique genome IDs
        genome_ids = set()
        genome_files = defaultdict(list)
        
        for bf in blast_files:
            basename = os.path.basename(bf)
            genome_id = basename.replace('_genes_blast.txt', '')
            genome_ids.add(genome_id)
            genome_files[genome_id].append(bf)
            
            # Also check for noncoding file
            noncoding_file = bf.replace('_genes_blast.txt', '_noncoding_blast.txt')
            if os.path.exists(noncoding_file):
                genome_files[genome_id].append(noncoding_file)
        
        print(f"Found {len(genome_ids)} genomes with BLAST results")
        print("Analyzing operon structure...")
        
        # Process each genome
        results = []
        operon_types = Counter()
        
        for i, genome_id in enumerate(sorted(genome_ids)):
            if i % 500 == 0:
                print(f"  Processed {i}/{len(genome_ids)} genomes...")
            
            result = self.analyze_genome_operon(genome_id, genome_files[genome_id])
            results.append(result)
            operon_types[result['operon_type']] += 1
        
        # Create DataFrame
        self.results_df = pd.DataFrame(results)
        
        # Save results
        output_file = os.path.join(self.output_dir, 'operon_order_analysis.csv')
        self.results_df.to_csv(output_file, index=False)
        print(f"\nDetailed results saved to: {output_file}")
        
        return True
    
    def generate_summary_statistics(self):
        """Generate and save summary statistics."""
        
        if self.results_df is None:
            print("Error: No analysis results available")
            return
        
        df = self.results_df
        total_genomes = len(df)
        
        print("\n" + "="*60)
        print("Summary Statistics")
        print("="*60)
        
        # Basic statistics
        complete_operons = df['complete_operon'].sum()
        print(f"Total genomes analyzed: {total_genomes}")
        print(f"Complete operons (7 genes): {complete_operons} ({100*complete_operons/total_genomes:.1f}%)")
        print()
        
        # Gene presence statistics
        print("Gene presence across genomes:")
        gene_counts = Counter()
        for genes in df['genes_found']:
            for gene in genes:
                gene_counts[gene] += 1
        
        for gene in CANONICAL_ORDER:
            count = gene_counts.get(gene, 0)
            print(f"  {gene} ({GENE_FUNCTIONS[gene]}): {count} ({100*count/total_genomes:.1f}%)")
        
        print()
        
        # Operon organization statistics
        print("Operon organization:")
        same_contig = df['same_contig'].sum()
        same_strand = df['same_strand'].sum()
        canonical = df['canonical_order'].sum()
        
        print(f"  On same contig: {same_contig} ({100*same_contig/total_genomes:.1f}%)")
        print(f"  On same strand: {same_strand} ({100*same_strand/total_genomes:.1f}%)")
        print(f"  Canonical gene order: {canonical} ({100*canonical/total_genomes:.1f}%)")
        print(f"  With inversions: {df['has_inversion'].sum()} ({100*df['has_inversion'].sum()/total_genomes:.1f}%)")
        print(f"  With rearrangements: {df['has_rearrangement'].sum()} ({100*df['has_rearrangement'].sum()/total_genomes:.1f}%)")
        
        # Save summary
        summary_file = os.path.join(self.output_dir, 'operon_order_summary.txt')
        with open(summary_file, 'w') as f:
            f.write("Operon Order and Synteny Analysis Summary\n")
            f.write("="*60 + "\n\n")
            f.write(f"Total genomes analyzed: {total_genomes}\n")
            f.write(f"Complete operons (7 genes): {complete_operons} ({100*complete_operons/total_genomes:.1f}%)\n\n")
            
            f.write("Gene presence:\n")
            for gene in CANONICAL_ORDER:
                count = gene_counts.get(gene, 0)
                f.write(f"  {gene}: {count} ({100*count/total_genomes:.1f}%)\n")
            
            f.write(f"\nOperon organization:\n")
            f.write(f"  Same contig: {same_contig} ({100*same_contig/total_genomes:.1f}%)\n")
            f.write(f"  Same strand: {same_strand} ({100*same_strand/total_genomes:.1f}%)\n")
            f.write(f"  Canonical order: {canonical} ({100*canonical/total_genomes:.1f}%)\n")
        
        print(f"\nSummary saved to: {summary_file}")
        
        # Create gene co-occurrence matrix
        print("\nCalculating gene co-occurrence...")
        cooccurrence = pd.DataFrame(0, index=CANONICAL_ORDER, columns=CANONICAL_ORDER)
        
        for genes in df['genes_found']:
            for gene1 in genes:
                for gene2 in genes:
                    if gene1 in CANONICAL_ORDER and gene2 in CANONICAL_ORDER:
                        cooccurrence.loc[gene1, gene2] += 1
        
        cooccurrence_file = os.path.join(self.output_dir, 'gene_cooccurrence_matrix.csv')
        cooccurrence.to_csv(cooccurrence_file)
        print(f"Gene co-occurrence matrix saved to: {cooccurrence_file}")
    
    def identify_non_canonical_orders(self):
        """Identify and report non-canonical gene orders."""
        
        if self.results_df is None:
            print("Error: No analysis results available")
            return
        
        df = self.results_df
        
        # Filter for complete operons
        complete_operons = df[df['complete_operon'] == True]
        print(f"\nFound {len(complete_operons)} genomes with complete operons")
        
        # Find non-canonical orders
        non_canonical = []
        canonical_count = 0
        
        for idx, row in complete_operons.iterrows():
            genome_id = row['genome_id']
            gene_order = row['gene_order']
            
            # Parse gene order if stored as string
            if isinstance(gene_order, str):
                try:
                    gene_order = eval(gene_order)
                except:
                    continue
            
            # Check if it matches canonical order
            if gene_order == CANONICAL_ORDER or gene_order == CANONICAL_ORDER[::-1]:
                canonical_count += 1
            else:
                non_canonical.append({
                    'genome_id': genome_id,
                    'gene_order': gene_order,
                    'same_contig': row['same_contig'],
                    'same_strand': row['same_strand'],
                    'synteny_score': row['synteny_score']
                })
        
        print(f"  Canonical order: {canonical_count}")
        print(f"  Non-canonical order: {len(non_canonical)}")
        
        # Save non-canonical orders
        output_file = os.path.join(self.output_dir, 'non_canonical_gene_orders.txt')
        with open(output_file, 'w') as f:
            f.write("Non-canonical Operon Gene Orders\n")
            f.write("="*60 + "\n\n")
            f.write(f"Total complete operons: {len(complete_operons)}\n")
            f.write(f"Canonical order: {canonical_count}\n")
            f.write(f"Non-canonical order: {len(non_canonical)}\n\n")
            
            if non_canonical:
                f.write("Detailed non-canonical orders:\n")
                f.write("-"*40 + "\n")
                for item in non_canonical:
                    f.write(f"\nGenome: {item['genome_id']}\n")
                    f.write(f"Gene order: {' -> '.join(item['gene_order'])}\n")
                    f.write(f"Same contig: {item['same_contig']}\n")
                    f.write(f"Same strand: {item['same_strand']}\n")
                    f.write(f"Synteny score: {item['synteny_score']:.2f}\n")
        
        print(f"Non-canonical orders saved to: {output_file}")
    
    def create_visualizations(self):
        """Create all visualizations."""
        
        if self.results_df is None:
            print("Error: No analysis results available")
            return
        
        df = self.results_df
        
        # 1. Operon type distribution pie chart
        self._create_operon_type_chart(df)
        
        # 2. Gene order visualization
        self._create_gene_order_pdf(df)
        
        # 3. Synteny heatmap
        self._create_synteny_heatmap(df)
    
    def _create_operon_type_chart(self, df):
        """Create horizontal bar chart of operon types for better readability."""
        
        operon_types = df['operon_type'].value_counts()
        total_count = len(df)
        
        # Create horizontal bar chart
        fig, ax = plt.subplots(figsize=(12, 6))
        
        # Sort by count
        operon_types = operon_types.sort_values(ascending=True)
        
        # Define colors for each type
        color_map = {
            'Complete and conserved': '#2E7D32',
            'Complete with rearrangement': '#FFA726',
            'Complete but fragmented': '#FF7043',
            'Near-complete': '#42A5F5',
            'Partial': '#BDBDBD',
            'No operon': '#F5F5F5'
        }
        colors = [color_map.get(x, '#999999') for x in operon_types.index]
        
        # Create bars
        y_pos = np.arange(len(operon_types))
        bars = ax.barh(y_pos, operon_types.values, color=colors, 
                      alpha=0.8, edgecolor='black', linewidth=1)
        
        # Add value labels
        for i, (bar, count) in enumerate(zip(bars, operon_types.values)):
            percentage = count / total_count * 100
            ax.text(bar.get_width() + max(operon_types.values)*0.01, 
                   bar.get_y() + bar.get_height()/2,
                   f'{count:,} ({percentage:.1f}%)', 
                   va='center', fontsize=10, fontweight='bold')
        
        # Customize plot
        ax.set_yticks(y_pos)
        ax.set_yticklabels(operon_types.index, fontsize=11)
        ax.set_xlabel('Number of Genomes', fontsize=12)
        ax.set_title('Distribution of Operon Organization Types', 
                    fontsize=14, fontweight='bold')
        ax.set_xlim(0, max(operon_types.values) * 1.15)
        
        # Add grid
        ax.grid(axis='x', alpha=0.3, linestyle='--')
        ax.set_axisbelow(True)
        
        plt.tight_layout()
        
        # Save figure
        output_file = os.path.join(self.output_dir, 'operon_organization_types.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Operon type chart saved to: {output_file}")
    
    def _create_gene_order_pdf(self, df):
        """Create PDF visualization of gene orders."""
        
        # Select representative examples
        examples = []
        operon_types = ['Complete and conserved', 'Complete with rearrangement', 
                       'Complete but fragmented', 'Near-complete', 'Partial']
        
        for otype in operon_types:
            type_df = df[df['operon_type'] == otype]
            if len(type_df) > 0:
                examples.append(type_df.iloc[0])
        
        if not examples:
            print("No examples found for gene order visualization")
            return
        
        # Create figure
        fig, axes = plt.subplots(len(examples) + 1, 1, figsize=(12, 2 * (len(examples) + 1)))
        
        # Plot canonical order first
        ax = axes[0] if len(examples) > 0 else axes
        self._plot_gene_order(ax, CANONICAL_ORDER, "Canonical Order", all_forward=True)
        
        # Plot examples
        for i, example in enumerate(examples):
            ax = axes[i + 1]
            gene_order = example['gene_order']
            if isinstance(gene_order, str):
                try:
                    gene_order = eval(gene_order)
                except:
                    continue
            
            # Check if fragmented
            is_fragmented = example.get('single_contig', True) == False
            missing_genes = list(set(CANONICAL_ORDER) - set(gene_order)) if len(gene_order) < 7 else []
            
            title = f"{example['operon_type']} ({example['genome_id'][:15]}...)"
            self._plot_gene_order(ax, gene_order, title, 
                                all_forward=example.get('same_strand', False),
                                genome_id=None,  # Already in title
                                is_fragmented=is_fragmented,
                                missing_genes=missing_genes)
        
        # Add overall title with explanation
        plt.suptitle('Operon Gene Order Patterns\n'
                    '"Conserved" = maintains canonical gene order | "//" = genes on different contigs',
                    fontsize=12, fontweight='bold', y=0.98)
        
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Leave space for suptitle
        
        # Save figure
        output_file = os.path.join(self.output_dir, 'gene_order_visualization.pdf')
        plt.savefig(output_file, format='pdf', bbox_inches='tight')
        plt.close()
        
        print(f"Gene order visualization saved to: {output_file}")
    
    def _plot_gene_order(self, ax, gene_order, title, all_forward=True,
                        genome_id=None, is_fragmented=False, missing_genes=None):
        """Plot a single gene order with fragmentation and missing gene indicators."""
        
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 1)
        ax.set_title(title, fontsize=10)
        ax.axis('off')
        
        # Calculate positions
        n_genes = len(gene_order)
        if n_genes == 0:
            return
        
        # Use fixed width for consistent visualization
        width = 7.0 / 7  # Always use width for 7 genes
        start_x = 1.0
        
        for i, gene in enumerate(gene_order):
            x = start_x + i * width
            color = GENE_COLORS.get(gene, '#CCCCCC')
            
            # Draw gene box
            rect = FancyBboxPatch((x, 0.3), width * 0.85, 0.4,
                                  boxstyle="round,pad=0.02",
                                  linewidth=1, 
                                  edgecolor='black',
                                  facecolor=color,
                                  alpha=0.9)
            ax.add_patch(rect)
            
            # Add gene name
            ax.text(x + width * 0.425, 0.5, gene,
                   ha='center', va='center',
                   fontsize=8, fontweight='bold')
            
            # Add fragmentation indicator in the middle
            if is_fragmented and i == n_genes // 2:
                ax.text(x - width * 0.2, 0.5, '//', fontsize=11, 
                       ha='center', va='center',
                       fontweight='bold', color='red')
            
            # Add arrow between genes
            if all_forward and i < n_genes - 1:
                ax.arrow(x + width * 0.8, 0.5, width * 0.08, 0,
                        head_width=0.04, head_length=0.02,
                        fc='black', ec='black', alpha=0.7)
        
        # Show missing genes for partial operons
        if missing_genes and len(missing_genes) > 0:
            x_end = start_x + n_genes * width
            ax.text(x_end + 0.2, 0.5, f'... ({len(missing_genes)} missing)', 
                   fontsize=9, ha='left', va='center', 
                   style='italic', color='gray')
    
    def _create_synteny_heatmap(self, df):
        """Skip creating the co-occurrence heatmap since all values are ~1.0.
        
        The co-occurrence matrix shows that when the operon is present,
        all genes co-occur nearly 100% of the time, making the heatmap
        uninformative (all cells show ~1.0). The CSV matrix is still
        generated in analyze_operon_order() for reference.
        """
        # Heatmap creation removed - not informative when all values are ~1.0
        pass
    
    def generate_manuscript_statistics(self):
        """Generate comprehensive manuscript statistics."""
        
        if self.results_df is None:
            print("Error: No analysis results available")
            return
        
        df = self.results_df
        total_genomes = len(df)
        
        output_lines = []
        output_lines.append("="*60)
        output_lines.append("OPERON ORDER ANALYSIS - MANUSCRIPT STATISTICS")
        output_lines.append("="*60)
        output_lines.append("")
        
        # Dataset overview
        output_lines.append("Dataset Overview:")
        output_lines.append(f"  Total genomes analyzed: {total_genomes:,}")
        output_lines.append("")
        
        # Operon completeness
        complete_operons = df['complete_operon'].sum()
        output_lines.append("Operon Completeness:")
        output_lines.append(f"  Complete operons (all 7 genes): {complete_operons:,} ({100*complete_operons/total_genomes:.1f}%)")
        
        # Gene count distribution
        gene_counts = df['n_genes_found'].value_counts().sort_index(ascending=False)
        output_lines.append(f"  Genomes by gene count:")
        for n_genes, count in gene_counts.items():
            output_lines.append(f"    {n_genes} genes: {count:,} ({100*count/total_genomes:.1f}%)")
        output_lines.append("")
        
        # Individual gene presence
        output_lines.append("Individual Gene Presence:")
        gene_presence = Counter()
        for genes in df['genes_found']:
            for gene in genes:
                gene_presence[gene] += 1
        
        for gene in CANONICAL_ORDER:
            count = gene_presence.get(gene, 0)
            pct = 100 * count / total_genomes
            output_lines.append(f"  {gene}: {count:,} genomes ({pct:.1f}%)")
        output_lines.append("")
        
        # Synteny conservation
        output_lines.append("Synteny Conservation (for complete operons):")
        complete_df = df[df['complete_operon'] == True]
        if len(complete_df) > 0:
            canonical = complete_df['canonical_order'].sum()
            same_contig = complete_df['same_contig'].sum()
            same_strand = complete_df['same_strand'].sum()
            
            output_lines.append(f"  Canonical gene order: {canonical:,} ({100*canonical/len(complete_df):.1f}%)")
            output_lines.append(f"  All genes on same contig: {same_contig:,} ({100*same_contig/len(complete_df):.1f}%)")
            output_lines.append(f"  All genes on same strand: {same_strand:,} ({100*same_strand/len(complete_df):.1f}%)")
            
            # Average synteny score
            mean_synteny = complete_df['synteny_score'].mean()
            output_lines.append(f"  Mean synteny score: {mean_synteny:.3f}")
        output_lines.append("")
        
        # Operon organization types
        output_lines.append("Operon Organization Types:")
        type_counts = df['operon_type'].value_counts()
        for otype, count in type_counts.items():
            output_lines.append(f"  {otype}: {count:,} ({100*count/total_genomes:.1f}%)")
        output_lines.append("")
        
        # Gene distance statistics
        valid_distances = df['mean_gene_distance'].dropna()
        if len(valid_distances) > 0:
            output_lines.append("Inter-gene Distance Statistics:")
            output_lines.append(f"  Mean distance: {valid_distances.mean():.0f} bp")
            output_lines.append(f"  Median distance: {valid_distances.median():.0f} bp")
            output_lines.append(f"  Min distance: {valid_distances.min():.0f} bp")
            output_lines.append(f"  Max distance: {df['max_gene_distance'].max():.0f} bp")
            output_lines.append("")
        
        # Rearrangements and inversions
        output_lines.append("Structural Variations:")
        inversions = df['has_inversion'].sum()
        rearrangements = df['has_rearrangement'].sum()
        output_lines.append(f"  Operons with inversions: {inversions:,} ({100*inversions/total_genomes:.1f}%)")
        output_lines.append(f"  Operons with rearrangements: {rearrangements:,} ({100*rearrangements/total_genomes:.1f}%)")
        output_lines.append("")
        
        # Save to file
        output_file = os.path.join(self.output_dir, 'manuscript_stats.txt')
        with open(output_file, 'w') as f:
            f.write('\n'.join(output_lines))
        
        # Also print to console
        print('\n'.join(output_lines))
        print(f"\nManuscript statistics saved to: {output_file}")


def stratify_by_metadata(order_summary_file, metadata_file, output_dir):
    """Stratify operon order results by metadata categories."""
    
    if not Path(metadata_file).exists():
        print(f"Metadata file not found: {metadata_file}, skipping stratification")
        return
    
    if not Path(order_summary_file).exists():
        print(f"Order summary file not found: {order_summary_file}, skipping stratification")
        return
    
    print("\n" + "="*60)
    print("Performing metadata stratification analysis")
    print("="*60)
    
    # Load operon order results
    order_df = pd.read_csv(order_summary_file)
    
    # Extract barcode from genome_id
    order_df['barcode'] = order_df['genome_id'].str.replace('_AS.result', '').str.replace('.result', '')
    
    # Load metadata
    meta_df = pd.read_csv(metadata_file, sep='\t')
    
    # Create mapping from both Barcode and AS_barcode
    meta_mapping = {}
    for _, row in meta_df.iterrows():
        barcode = row['Barcode']
        meta_mapping[barcode] = row
        if 'AS_barcode' in row and pd.notna(row['AS_barcode']) and row['AS_barcode'] != 'ND':
            as_barcode = row['AS_barcode'].replace('_AS', '')
            meta_mapping[as_barcode] = row
    
    # Merge metadata
    merged_data = []
    unmapped = []
    for _, row in order_df.iterrows():
        barcode = row['barcode']
        if barcode in meta_mapping:
            meta_row = meta_mapping[barcode]
            merged_row = row.to_dict()
            merged_row['Source Niche'] = meta_row.get('Source Niche', 'Unknown')
            merged_row['Source Details'] = meta_row.get('Source Details', 'Unknown')
            merged_row['Country'] = meta_row.get('Country', 'Unknown')
            merged_data.append(merged_row)
        else:
            unmapped.append(barcode)
    
    if not merged_data:
        print("Warning: No genomes could be matched to metadata")
        return
    
    merged_df = pd.DataFrame(merged_data)
    print(f"Matched {len(merged_df)} genomes to metadata ({len(unmapped)} unmapped)")
    
    # Analyze patterns by metadata categories
    output_lines = []
    output_lines.append("="*60)
    output_lines.append("OPERON ORDER STRATIFIED BY METADATA")
    output_lines.append("="*60)
    output_lines.append("")
    
    # 1. By Source Niche
    output_lines.append("OPERON COMPLETENESS BY SOURCE NICHE:")
    output_lines.append("-"*40)
    
    niche_results = []
    for niche in sorted(merged_df['Source Niche'].unique()):
        if pd.isna(niche) or niche == 'ND':
            continue
        niche_df = merged_df[merged_df['Source Niche'] == niche]
        total = len(niche_df)
        
        # Count complete operons (all 7 genes)
        complete = (niche_df['num_genes'] == 7).sum()
        complete_pct = 100 * complete / total if total > 0 else 0
        
        # Count canonical order
        canonical = (niche_df['is_canonical'] == True).sum() if 'is_canonical' in niche_df.columns else 0
        canonical_pct = 100 * canonical / total if total > 0 else 0
        
        # Average number of genes
        avg_genes = niche_df['num_genes'].mean()
        
        niche_results.append({
            'Source Niche': niche,
            'Total Genomes': total,
            'Complete Operons': complete,
            'Complete %': round(complete_pct, 1),
            'Canonical Order': canonical,
            'Canonical %': round(canonical_pct, 1),
            'Avg Genes': round(avg_genes, 2)
        })
        
        output_lines.append(f"  {niche}: n={total}")
        output_lines.append(f"    Complete operons: {complete} ({complete_pct:.1f}%)")
        output_lines.append(f"    Canonical order: {canonical} ({canonical_pct:.1f}%)")
        output_lines.append(f"    Average genes: {avg_genes:.2f}")
    
    if niche_results:
        niche_df_out = pd.DataFrame(niche_results)
        niche_df_out.to_csv(os.path.join(output_dir, "operon_order_by_source_niche.tsv"), sep='\t', index=False)
    
    # 2. By Country (Top 10)
    output_lines.append("")
    output_lines.append("OPERON COMPLETENESS BY COUNTRY (Top 10):")
    output_lines.append("-"*40)
    
    country_results = []
    country_counts = merged_df['Country'].value_counts()
    for country in country_counts.head(10).index:
        if pd.isna(country) or country == 'ND' or country == 'Unknown':
            continue
        country_df = merged_df[merged_df['Country'] == country]
        total = len(country_df)
        
        complete = (country_df['num_genes'] == 7).sum()
        complete_pct = 100 * complete / total if total > 0 else 0
        canonical = (country_df['is_canonical'] == True).sum() if 'is_canonical' in country_df.columns else 0
        canonical_pct = 100 * canonical / total if total > 0 else 0
        
        country_results.append({
            'Country': country,
            'Total': total,
            'Complete': complete,
            'Complete%': round(complete_pct, 1),
            'Canonical': canonical,
            'Canonical%': round(canonical_pct, 1)
        })
        
        output_lines.append(f"  {country}: n={total}, Complete={complete_pct:.1f}%, Canonical={canonical_pct:.1f}%")
    
    if country_results:
        country_df_out = pd.DataFrame(country_results)
        country_df_out.to_csv(os.path.join(output_dir, "operon_order_by_country.tsv"), sep='\t', index=False)
    
    # 3. Gene presence patterns by Source Niche
    output_lines.append("")
    output_lines.append("GENE PRESENCE BY SOURCE NICHE:")
    output_lines.append("-"*40)
    
    gene_presence_results = []
    for niche in sorted(merged_df['Source Niche'].unique()):
        if pd.isna(niche) or niche == 'ND':
            continue
        niche_df = merged_df[merged_df['Source Niche'] == niche]
        
        gene_counts = {}
        for gene in CANONICAL_ORDER:
            if f'has_{gene}' in niche_df.columns:
                present = niche_df[f'has_{gene}'].sum()
                gene_counts[gene] = present
            elif 'genes_present' in niche_df.columns:
                # Parse genes_present string if stored that way
                present = niche_df['genes_present'].str.contains(gene).sum()
                gene_counts[gene] = present
        
        if gene_counts:
            result = {'Source Niche': niche, 'Total': len(niche_df)}
            result.update(gene_counts)
            gene_presence_results.append(result)
            
            output_lines.append(f"  {niche} (n={len(niche_df)}):")
            gene_str = ", ".join([f"{g}:{c}" for g, c in gene_counts.items()])
            output_lines.append(f"    {gene_str}")
    
    if gene_presence_results:
        gene_df = pd.DataFrame(gene_presence_results)
        gene_df.to_csv(os.path.join(output_dir, "gene_presence_by_source_niche.tsv"), sep='\t', index=False)
    
    # 4. Summary statistics
    output_lines.append("")
    output_lines.append("SUMMARY:")
    output_lines.append("-"*40)
    output_lines.append(f"Total genomes with metadata: {len(merged_df)}")
    output_lines.append(f"Unique source niches: {merged_df['Source Niche'].nunique()}")
    output_lines.append(f"Unique countries: {merged_df['Country'].nunique()}")
    
    overall_complete = (merged_df['num_genes'] == 7).sum()
    overall_complete_pct = 100 * overall_complete / len(merged_df)
    output_lines.append(f"Overall complete operons: {overall_complete} ({overall_complete_pct:.1f}%)")
    
    if 'is_canonical' in merged_df.columns:
        overall_canonical = merged_df['is_canonical'].sum()
        overall_canonical_pct = 100 * overall_canonical / len(merged_df)
        output_lines.append(f"Overall canonical order: {overall_canonical} ({overall_canonical_pct:.1f}%)")
    
    # Write summary to file
    summary_path = os.path.join(output_dir, "operon_order_metadata_stratification.txt")
    with open(summary_path, 'w') as f:
        f.write('\n'.join(output_lines))
    
    print('\n'.join(output_lines))
    print(f"\nStratification results saved to:")
    print(f"  - operon_order_by_source_niche.tsv")
    print(f"  - operon_order_by_country.tsv")
    print(f"  - gene_presence_by_source_niche.tsv")
    print(f"  - operon_order_metadata_stratification.txt")


def main():
    """Main function to run the complete pipeline."""
    
    parser = argparse.ArgumentParser(description="Unified operon order and synteny analysis pipeline")
    parser.add_argument("--blast-dir",
                       default="../03_blast_search/output/blast_results",
                       help="Directory containing BLAST results")
    parser.add_argument("--output-dir", 
                       default="output",
                       help="Output directory")
    parser.add_argument("--min-identity", 
                       type=float, 
                       default=80,
                       help="Minimum sequence identity threshold")
    parser.add_argument("--steps",
                       nargs='+',
                       choices=['analyze', 'visualize', 'stats', 'all'],
                       default=['all'],
                       help="Which analysis steps to run")
    parser.add_argument("--metadata",
                       default="../00_annotation/8587_Efs_metadata_ASbarcode.txt",
                       help="Metadata file for stratification analysis")
    
    args = parser.parse_args()
    
    # Initialize analyzer
    analyzer = OperonOrderAnalyzer(args.blast_dir, args.output_dir, args.min_identity)
    
    # Determine which steps to run
    if 'all' in args.steps:
        steps = ['analyze', 'visualize', 'stats']
    else:
        steps = args.steps
    
    # Run analysis
    if 'analyze' in steps:
        print("\n" + "="*60)
        print("Step 1: Running operon order analysis")
        print("="*60)
        
        if not analyzer.run_analysis():
            print("Error: Analysis failed")
            return 1
        
        analyzer.generate_summary_statistics()
        analyzer.identify_non_canonical_orders()
    
    # Load existing results if not analyzing
    elif 'visualize' in steps or 'stats' in steps:
        results_file = os.path.join(args.output_dir, 'operon_order_analysis.csv')
        if os.path.exists(results_file):
            analyzer.results_df = pd.read_csv(results_file)
        else:
            print(f"Error: Results file not found: {results_file}")
            print("Please run analysis first")
            return 1
    
    # Create visualizations
    if 'visualize' in steps:
        print("\n" + "="*60)
        print("Step 2: Creating visualizations")
        print("="*60)
        analyzer.create_visualizations()
    
    # Generate manuscript statistics
    if 'stats' in steps:
        print("\n" + "="*60)
        print("Step 3: Generating manuscript statistics")
        print("="*60)
        analyzer.generate_manuscript_statistics()
    
    # Perform metadata stratification if metadata file exists
    if args.metadata and Path(args.metadata).exists():
        order_summary_file = os.path.join(args.output_dir, "operon_order_summary.csv")
        if Path(order_summary_file).exists():
            stratify_by_metadata(order_summary_file, args.metadata, args.output_dir)
        else:
            print(f"\nSkipping metadata stratification - order summary not found: {order_summary_file}")
    
    print("\n" + "="*60)
    print("Pipeline complete!")
    print("="*60)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())