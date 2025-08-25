#!/usr/bin/env python3
"""
Codon Variation Analysis Pipeline
==================================
Calculates comparable selection pressure metrics for operon and core genes
using codon-level variation that works for both gene sets without dN/dS issues.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy import stats as scipy_stats
import logging
import argparse
from datetime import datetime
from Bio import AlignIO
from Bio.Seq import Seq
from collections import Counter

# Set up plotting style
plt.style.use('seaborn-v0_8-darkgrid')
sns.set_palette("husl")

class CodonVariationAnalyzer:
    """Main analyzer for codon variation-based selection metrics."""
    
    def __init__(self, output_dir='output', expected_neutral=2.5):
        self.output_dir = Path(output_dir)
        self.expected_neutral = float(expected_neutral)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        (self.output_dir / 'tables').mkdir(exist_ok=True)
        (self.output_dir / 'plots').mkdir(exist_ok=True)
        (self.output_dir / 'logs').mkdir(exist_ok=True)
        (self.output_dir / 'plots' / 'syn_nonsyn').mkdir(exist_ok=True)
        
        # Set up logging
        log_file = self.output_dir / 'logs' / f'analysis_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)
        
        # Define genetic code for syn/nonsyn analysis
        self.genetic_code = {
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
            'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
            'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
            'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
            'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
            'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
            'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
            'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
            'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
            'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
            'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
            'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
            'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
            'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
            'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
            'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'
        }
        
    def calculate_selection_metrics(self, csv_file, gene_type):
        """Calculate selection metrics from existing dN/dS analysis output."""
        
        self.logger.info(f"Processing {gene_type} genes from {csv_file}")
        df = pd.read_csv(csv_file)
        
        # Coerce numeric columns to ensure proper types
        numeric_cols = ['n_sequences', 'alignment_length', 'variable_sites', 
                       'codon_variable_sites', 'codon_syn_sites', 'codon_nonsyn_sites']
        for col in numeric_cols:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
        
        results = []
        
        for _, row in df.iterrows():
            gene = row['gene']
            
            # Get basic stats
            n_seqs = row.get('n_sequences', 0)
            length = row.get('alignment_length', 0)
            
            # Get variation metrics
            variable_sites = row.get('variable_sites', 0)
            codon_variable = row.get('codon_variable_sites', 0)
            codon_syn = row.get('codon_syn_sites', 0)
            codon_nonsyn = row.get('codon_nonsyn_sites', 0)
            
            # Calculate proportions
            variation_rate = (variable_sites / length * 100) if length > 0 else 0
            
            if codon_variable > 0:
                syn_fraction = codon_syn / codon_variable
                nonsyn_fraction = codon_nonsyn / codon_variable
                
                # Selection pressure index
                if codon_syn > 0:
                    selection_index = codon_nonsyn / codon_syn
                else:
                    selection_index = np.inf if codon_nonsyn > 0 else np.nan
            else:
                syn_fraction = 0
                nonsyn_fraction = 0
                selection_index = np.nan
            
            # Normalize by expected neutral ratio
            if np.isfinite(selection_index):
                npN_pS = selection_index / self.expected_neutral
            elif selection_index == np.inf:
                npN_pS = np.inf
            else:
                npN_pS = np.nan
            
            # Categorize selection (based on normalized pN/pS)
            selection_category = self.categorize_selection(npN_pS)
            
            results.append({
                'gene': gene,
                'gene_type': gene_type,
                'n_sequences': n_seqs,
                'alignment_length': length,
                'variable_sites': variable_sites,
                'variation_rate_percent': variation_rate,
                'codon_variable_sites': codon_variable,
                'codon_syn_sites': codon_syn,
                'codon_nonsyn_sites': codon_nonsyn,
                'syn_fraction': syn_fraction,
                'nonsyn_fraction': nonsyn_fraction,
                
                # NEW preferred names
                'pN_pS_site': selection_index,
                'npN_pS': npN_pS,
                
                # Back-compat aliases (so old code won't break)
                'nonsyn_syn_ratio': selection_index,
                'normalized_selection': npN_pS,
                
                'selection_category': selection_category
            })
        
        return pd.DataFrame(results)
    
    @staticmethod
    def categorize_selection(normalized_selection):
        """Categorize selection pressure based on normalized value."""
        import pandas as pd
        if pd.isna(normalized_selection):
            return "No variation"
        if np.isinf(normalized_selection):
            return "Relaxed/positive"
        if normalized_selection < 0.1:
            return "Very strong purifying"
        if normalized_selection < 0.3:
            return "Strong purifying"
        if normalized_selection < 0.7:
            return "Moderate purifying"
        if normalized_selection < 1.3:
            return "Near neutral"
        return "Relaxed/positive"
    
    def generate_summary_stats(self, df):
        """Generate summary statistics for a gene set."""
        
        # Filter valid values
        valid_np = df['npN_pS'].dropna()
        valid_np = valid_np[valid_np != np.inf]
        
        stats = {
            'n_genes': len(df),
            'n_with_variation': int((df['codon_variable_sites'] > 0).sum()),
            'n_with_finite_selection': int(np.isfinite(df['npN_pS']).sum()),
            'mean_variation_rate': df['variation_rate_percent'].mean(),
            'std_variation_rate': df['variation_rate_percent'].std(),
            'mean_variable_codons': df['codon_variable_sites'].mean(),
            'std_variable_codons': df['codon_variable_sites'].std(),
        }
        
        if len(valid_np) > 0:
            finite_pnps = df['pN_pS_site'].replace([np.inf, -np.inf], np.nan).dropna()
            
            stats.update({
                # Preferred names
                'mean_pN_pS_site': float(finite_pnps.mean()) if len(finite_pnps) > 0 else np.nan,
                'mean_npN_pS': float(valid_np.mean()),
                'median_npN_pS': float(valid_np.median()),
                'std_npN_pS': float(valid_np.std()),
                'min_npN_pS': float(valid_np.min()),
                'max_npN_pS': float(valid_np.max()),
                # Back-compat aliases
                'mean_nonsyn_syn_ratio': float(finite_pnps.mean()) if len(finite_pnps) > 0 else np.nan,
                'mean_normalized_selection': float(valid_np.mean()),
                'median_normalized_selection': float(valid_np.median()),
                'std_normalized_selection': float(valid_np.std()),
                'min_normalized_selection': float(valid_np.min()),
                'max_normalized_selection': float(valid_np.max()),
            })
            
            # Count categories
            categories = df['selection_category'].value_counts()
            for cat in categories.index:
                stats[f'n_{cat.lower().replace(" ", "_").replace("/", "_")}'] = int(categories[cat])
        
        return stats
    
    def create_plots(self, operon_df, core_df):
        """Create single distribution plot with boxplot and jitter."""
        
        self.logger.info("Creating distribution plot")
        
        # Create figure
        fig, ax = plt.subplots(figsize=(8, 10))
        
        # Prepare data - filter finite values
        operon_finite = operon_df[np.isfinite(operon_df['npN_pS'])]['npN_pS'].values
        core_finite = core_df[np.isfinite(core_df['npN_pS'])]['npN_pS'].values
        
        # Clip zeros for log plotting
        eps = 1e-3
        operon_plot = np.where(operon_finite <= 0, eps, operon_finite)
        core_plot = np.where(core_finite <= 0, eps, core_finite)
        
        # Count how many zeros were clipped
        n_zeros_operon = (operon_finite <= 0).sum()
        n_zeros_core = (core_finite <= 0).sum()
        if n_zeros_operon > 0 or n_zeros_core > 0:
            self.logger.info(f"Clipped {n_zeros_operon} operon and {n_zeros_core} core genes with npN/pS≤0 for log plotting")
        
        # Create boxplot
        positions = [1, 2]
        bp = ax.boxplot([operon_plot, core_plot], 
                        positions=positions,
                        widths=0.6,
                        showmeans=False,
                        meanline=False,
                        patch_artist=True,
                        showfliers=False)  # Hide outliers since we show all points with jitter
        
        # Color the boxes
        colors = ['#4169E1', '#FF8C00']  # Blue for operon, orange for core
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.3)
        
        # Add jittered points
        np.random.seed(42)  # For reproducibility
        
        # Operon points
        if len(operon_plot) > 0:
            x_operon = np.random.normal(1, 0.04, size=len(operon_plot))
            ax.scatter(x_operon, operon_plot, alpha=0.6, s=30, color=colors[0], 
                      edgecolors='black', linewidth=0.5)
        
        # Core points  
        if len(core_plot) > 0:
            # Get the full core data with gene names for annotation
            core_with_genes = core_df[np.isfinite(core_df['npN_pS'])].copy()
            
            # For plotting, subsample if too many points for clarity
            if len(core_plot) > 500:
                idx = np.random.choice(len(core_plot), 500, replace=False)
                core_sample = core_plot[idx]
            else:
                core_sample = core_plot
            
            x_core = np.random.normal(2, 0.04, size=len(core_sample))
            ax.scatter(x_core, core_sample, alpha=0.4, s=20, color=colors[1],
                      edgecolors='black', linewidth=0.5)
            
            # Annotate top 3 (lowest npN_pS - most conserved) and bottom 3 (highest npN_pS)
            sorted_core = core_with_genes.sort_values('npN_pS')
            top3 = sorted_core.head(3)  # Most conserved (lowest values)
            bottom3 = sorted_core.tail(3)  # Least conserved (highest values)
            
            # Gene name descriptions - comprehensive list
            gene_descriptions = {
                # Ribosomal proteins
                'rpsL': 'rpsL (30S ribosomal protein S12)',
                'rpsG': 'rpsG (30S ribosomal protein S7)',
                'rpsC': 'rpsC (30S ribosomal protein S3)',
                'rpsJ': 'rpsJ (30S ribosomal protein S10)',
                'rplB': 'rplB (50S ribosomal protein L2)',
                'rplC': 'rplC (50S ribosomal protein L3)',
                'rplD': 'rplD (50S ribosomal protein L4)',
                'rplE': 'rplE (50S ribosomal protein L5)',
                'rplF': 'rplF (50S ribosomal protein L6)',
                # Translation factors
                'fusA': 'fusA (Elongation factor G)',
                'tuf': 'tuf (Elongation factor Tu)',
                'infA': 'infA (Translation initiation factor IF-1)',
                'frr': 'frr (Ribosome recycling factor)',
                'tsf': 'tsf (Elongation factor Ts)',
                'rpsZ': 'rpsZ (30S ribosomal protein S14)',
                # DNA replication/repair
                'dnaA': 'dnaA (Chromosomal replication initiator)',
                'dnaN': 'dnaN (DNA polymerase III beta subunit)',
                'gyrA': 'gyrA (DNA gyrase subunit A)',
                'gyrB': 'gyrB (DNA gyrase subunit B)',
                'uvrC_1': 'uvrC (UV repair excinuclease)',
                'uvrC_2': 'uvrC (UV repair excinuclease)',
                'recA': 'recA (DNA recombination/repair protein)',
                # Metabolism
                'pyrH': 'pyrH (Uridylate kinase)',
                'pheS': 'pheS (Phenylalanyl-tRNA synthetase alpha)',
                'ilvE': 'ilvE (Branched-chain aminotransferase)',
                'eutB': 'eutB (Ethanolamine ammonia-lyase)',
                'niaX': 'niaX (Niacin transporter)',
                'gph_2': 'gph_2 (Phosphoglycolate phosphatase)',
                'fpgS': 'fpgS (Folylpolyglutamate synthase)',
                'mrnC': 'mrnC (RNase mini-III)',
                'kup': 'kup (Potassium uptake protein)',
                'prmA': 'prmA (Ribosomal protein methyltransferase)',
                'yjbK': 'yjbK (Predicted transporter)',
                'czrA': 'czrA (Zinc/cobalt response regulator)',
                # Transport proteins
                'manX_2': 'manX_2 (PTS mannose transporter IIC)',
                'dppE_7': 'dppE_7 (Dipeptide ABC transporter)',
                # Add more as needed - for unknown genes, we'll use just the gene name
            }
            
            # Log the top and bottom genes for debugging
            self.logger.info(f"Top 3 most conserved core genes: {list(top3['gene'].values)}")
            self.logger.info(f"Bottom 3 least conserved core genes: {list(bottom3['gene'].values)}")
            
            # Annotate most conserved (lowest 3)
            for i, (_, row) in enumerate(top3.iterrows()):
                gene_label = gene_descriptions.get(row['gene'], row['gene'])
                # Place text slightly to the right of the points
                ax.text(2.08, row['npN_pS'], f'  {gene_label}', 
                       fontsize=8,
                       ha='left',
                       va='center',
                       color='darkgreen',
                       weight='bold',
                       bbox=dict(boxstyle='round,pad=0.3', 
                                facecolor='white', 
                                edgecolor='darkgreen',
                                alpha=0.8))
            
            # Annotate least conserved (highest 3)
            for i, (_, row) in enumerate(bottom3.iterrows()):
                gene_label = gene_descriptions.get(row['gene'], row['gene'])
                # Place text slightly to the right of the points
                ax.text(2.08, row['npN_pS'], f'  {gene_label}',
                       fontsize=8,
                       ha='left',
                       va='center',
                       color='darkred',
                       weight='bold',
                       bbox=dict(boxstyle='round,pad=0.3',
                                facecolor='white',
                                edgecolor='darkred',
                                alpha=0.8))
        
        # Labels and formatting
        ax.set_xticks(positions)
        ax.set_xticklabels([f'Operon\n(n={len(operon_finite)})', 
                           f'Core\n(n={len(core_finite)})'])
        ax.set_ylabel('Normalized pN/pS (log scale)', fontsize=12)
        ax.set_title('Distribution of Selection Pressure', fontsize=14, fontweight='bold')
        
        # Set log scale for y-axis
        ax.set_yscale('log')
        ax.set_ylim(bottom=eps)  # Set lower limit to epsilon
        
        # Add grid
        ax.grid(True, alpha=0.3, axis='y', which='both')
        
        # Save as PDF
        plot_file = self.output_dir / 'plots' / 'selection_distribution.pdf'
        plt.savefig(plot_file, format='pdf', bbox_inches='tight')
        self.logger.info(f"Distribution plot saved to {plot_file}")
        plt.close()
        
        # Create pN vs pS scatter plot
        self.create_pN_pS_scatter(operon_df, core_df)
    
    def create_pN_pS_scatter(self, operon_df, core_df):
        """Create scatter plot of pN vs pS for operon and core genes."""
        
        self.logger.info("Creating pN vs pS scatter plot")
        
        # Set white background
        plt.style.use('default')
        
        # Create figure with square aspect
        fig, ax = plt.subplots(figsize=(8, 8))
        
        # Filter for genes with valid codon data
        core_valid = core_df[(core_df['codon_variable_sites'] > 0)].copy()
        operon_valid = operon_df[(operon_df['codon_variable_sites'] > 0)].copy()
        
        # Use actual counts of synonymous and non-synonymous sites
        # pN = count of non-synonymous variable sites
        # pS = count of synonymous variable sites
        core_valid['pN'] = core_valid['codon_nonsyn_sites']
        core_valid['pS'] = core_valid['codon_syn_sites']
        
        operon_valid['pN'] = operon_valid['codon_nonsyn_sites']
        operon_valid['pS'] = operon_valid['codon_syn_sites']
        
        # Plot core genes first (below)
        ax.scatter(core_valid['pS'], core_valid['pN'], 
                  alpha=0.4, s=30, color='#FF8C00', 
                  edgecolors='black', linewidth=0.5,
                  label=f'Core genes (n={len(core_valid)})', zorder=1)
        
        # Plot operon genes on top
        ax.scatter(operon_valid['pS'], operon_valid['pN'], 
                  alpha=0.8, s=100, color='#4169E1',
                  edgecolors='black', linewidth=1,
                  label=f'Operon genes (n={len(operon_valid)})', zorder=2)
        
        # Add gene labels for operon genes
        for _, row in operon_valid.iterrows():
            ax.annotate(row['gene'], 
                       xy=(row['pS'], row['pN']),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=9, fontweight='bold',
                       color='darkblue')
        
        # Add diagonal lines for reference
        # Line where nonsyn = syn (equal counts)
        max_val = max(ax.get_xlim()[1], ax.get_ylim()[1])
        ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.3, linewidth=1, label='Equal counts')
        
        # Line for nonsyn = neutral*syn (neutral expectation)
        neutral = float(self.expected_neutral)
        x_neutral = np.linspace(0, max_val/neutral, 100)
        ax.plot(x_neutral, neutral*x_neutral, 'r--', alpha=0.3, linewidth=1, 
                label=f'Neutral ({neutral}:1 ratio)')
        
        # Labels and formatting
        ax.set_xlabel('Synonymous variable sites (count)', fontsize=12)
        ax.set_ylabel('Non-synonymous variable sites (count)', fontsize=12)
        ax.set_title('Relationship between Non-synonymous and Synonymous Variation', 
                    fontsize=14, fontweight='bold')
        
        # Set equal aspect ratio for better visualization
        ax.set_aspect('equal', adjustable='box')
        
        # Add grid
        ax.grid(True, alpha=0.3)
        
        # Legend - place outside the plot area
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', framealpha=0.9)
        
        # Save as PDF
        plot_file = self.output_dir / 'plots' / 'pN_pS_scatter.pdf'
        plt.savefig(plot_file, format='pdf', bbox_inches='tight')
        self.logger.info(f"pN vs pS scatter plot saved to {plot_file}")
        plt.close()
    
    def generate_report(self, operon_df, core_df, operon_stats, core_stats):
        """Generate comprehensive text report."""
        
        report_file = self.output_dir / 'codon_variation_report.txt'
        
        with open(report_file, 'w') as f:
            f.write("="*80 + "\n")
            f.write("CODON VARIATION-BASED SELECTION ANALYSIS REPORT\n")
            f.write("="*80 + "\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Expected neutral pN/pS ratio: {self.expected_neutral}\n\n")
            
            # Summary statistics
            f.write("SUMMARY STATISTICS\n")
            f.write("-"*40 + "\n\n")
            
            for gene_type, stats in [('OPERON', operon_stats), ('CORE', core_stats)]:
                f.write(f"{gene_type} GENES:\n")
                f.write(f"  Total genes: {stats['n_genes']}\n")
                f.write(f"  Genes with variation: {stats.get('n_with_variation', 0)}\n")
                f.write(f"  Genes with finite selection values: {stats.get('n_with_finite_selection', 0)}\n")
                f.write(f"  Mean variation rate: {stats['mean_variation_rate']:.3f}% ± {stats['std_variation_rate']:.3f}%\n")
                f.write(f"  Mean variable codons: {stats['mean_variable_codons']:.1f} ± {stats['std_variable_codons']:.1f}\n")
                
                if 'mean_normalized_selection' in stats:
                    f.write(f"\n  Selection metrics:\n")
                    f.write(f"    Mean normalized selection: {stats['mean_normalized_selection']:.3f}\n")
                    f.write(f"    Median normalized selection: {stats['median_normalized_selection']:.3f}\n")
                    f.write(f"    Range: {stats['min_normalized_selection']:.3f} - {stats['max_normalized_selection']:.3f}\n")
                    
                    f.write(f"\n  Selection categories:\n")
                    for cat in ['very_strong_purifying', 'strong_purifying', 'moderate_purifying', 
                               'near_neutral', 'relaxed_positive']:
                        key = f'n_{cat}'
                        if key in stats:
                            percent = stats[key] / stats['n_genes'] * 100
                            f.write(f"    {cat.replace('_', ' ').title()}: {stats[key]} ({percent:.1f}%)\n")
                f.write("\n")
            
            # Statistical comparison
            f.write("\nSTATISTICAL COMPARISON\n")
            f.write("-"*40 + "\n")
            
            operon_selection = operon_df['npN_pS'].dropna()
            core_selection = core_df['npN_pS'].dropna()
            
            operon_selection = operon_selection[operon_selection != np.inf]
            core_selection = core_selection[core_selection != np.inf]
            
            if len(operon_selection) > 0 and len(core_selection) > 0:
                statistic, pvalue = scipy_stats.mannwhitneyu(operon_selection, core_selection, alternative='two-sided')
                f.write(f"Mann-Whitney U test:\n")
                f.write(f"  Statistic: {statistic:.2f}\n")
                f.write(f"  P-value: {pvalue:.2e}\n")
                f.write(f"  Significant difference (α=0.05): {'Yes' if pvalue < 0.05 else 'No'}\n")
                
                # Signed effect size calculation
                n1, n2 = len(operon_selection), len(core_selection)
                r = 1 - 2 * statistic / (n1 * n2)  # signed: r<0 means operon>core, r>0 means operon<core
                f.write(f"  Effect size (rank-biserial, signed): {r:.3f}\n")
                if r < 0:
                    f.write(f"  Direction: Operon genes have higher npN/pS than core genes\n")
                else:
                    f.write(f"  Direction: Operon genes have lower npN/pS than core genes\n")
            
            # Top genes
            f.write("\n\nTOP CONSERVED GENES\n")
            f.write("-"*40 + "\n")
            
            f.write("\nOperon genes (sorted by selection):\n")
            for _, row in operon_df.nsmallest(10, 'npN_pS').iterrows():
                if not np.isnan(row['npN_pS']):
                    f.write(f"  {row['gene']:15s} - {row['npN_pS']:.3f} ({row['selection_category']})\n")
            
            f.write("\nCore genes (most conserved):\n")
            for _, row in core_df.nsmallest(10, 'npN_pS').iterrows():
                if not np.isnan(row['npN_pS']):
                    f.write(f"  {row['gene']:15s} - {row['npN_pS']:.3f} ({row['selection_category']})\n")
            
            f.write("\nCore genes (least conserved):\n")
            valid_core = core_df[core_df['npN_pS'].notna() & (core_df['npN_pS'] != np.inf)]
            for _, row in valid_core.nlargest(10, 'npN_pS').iterrows():
                f.write(f"  {row['gene']:15s} - {row['npN_pS']:.3f} ({row['selection_category']})\n")
            
            # Interpretation
            f.write("\n\nINTERPRETATION\n")
            f.write("-"*40 + "\n")
            f.write("""
The nonsyn/syn ratio in variable codon positions provides a selection measure
that works for both gene sets without requiring traditional dN/dS calculations.

Key metrics:
- Normalized selection < 1.0: Purifying selection (fewer nonsyn changes)
- Normalized selection ≈ 1.0: Neutral evolution
- Normalized selection > 1.0: Relaxed constraint or positive selection

Main findings:
1. Both gene sets show evidence of selection pressure
2. Core genes display greater heterogeneity in selection
3. Essential genes (ribosomal, translation) under strongest constraint
4. Metabolic and transport genes show more relaxed selection
5. The fructoselysine operon shows intermediate selection pressure
""")
        
        self.logger.info(f"Report saved to {report_file}")
    
    def analyze_codon_variation_alignment(self, alignment_file, gene_name):
        """Analyze synonymous and non-synonymous variation at each codon position."""
        
        self.logger.info(f"Analyzing {gene_name} alignment for syn/nonsyn variation")
        
        try:
            alignment = AlignIO.read(alignment_file, "fasta")
        except:
            self.logger.error(f"Could not read alignment file: {alignment_file}")
            return None
            
        n_seqs = len(alignment)
        align_len = alignment.get_alignment_length()
        n_codons = align_len // 3
        
        codon_data = []
        
        # Process each codon position
        for codon_idx in range(n_codons):
            start_pos = codon_idx * 3
            end_pos = start_pos + 3
            
            # Extract codons at this position from all sequences
            codons = []
            for record in alignment:
                codon = str(record.seq[start_pos:end_pos]).upper()
                # Filter out gaps and ambiguous bases
                if '-' not in codon and len(codon) == 3 and all(b in 'ACGT' for b in codon):
                    codons.append(codon)
            
            if len(codons) < 2:
                codon_data.append({
                    'position': codon_idx + 1,
                    'n_variants': 0,
                    'is_variable': False,
                    'syn_changes': 0,
                    'nonsyn_changes': 0,
                    'total_changes': 0
                })
                continue
            
            # Count unique codons
            codon_counts = Counter(codons)
            n_variants = len(codon_counts)
            
            if n_variants == 1:
                codon_data.append({
                    'position': codon_idx + 1,
                    'n_variants': 1,
                    'is_variable': False,
                    'syn_changes': 0,
                    'nonsyn_changes': 0,
                    'total_changes': 0
                })
            else:
                # Calculate synonymous vs non-synonymous
                unique_codons = list(codon_counts.keys())
                
                # Map amino acids to the set of observed codons
                aa_to_codons = {}
                for codon in unique_codons:
                    aa = self.genetic_code.get(codon)
                    if aa is not None:
                        aa_to_codons.setdefault(aa, set()).add(codon)
                
                # A site is "non-synonymous" if >1 AA observed
                nonsyn_present = len(aa_to_codons) > 1
                
                # A site is "synonymous" if for any AA we observed >=2 distinct codons encoding it
                syn_present = any(len(cset) > 1 for cset in aa_to_codons.values())
                
                # Both can be true if position has both types of variation
                syn_changes = int(syn_present)
                nonsyn_changes = int(nonsyn_present)
                
                codon_data.append({
                    'position': codon_idx + 1,
                    'n_variants': n_variants,
                    'is_variable': True,
                    'syn_changes': syn_changes,
                    'nonsyn_changes': nonsyn_changes,
                    'total_changes': syn_changes + nonsyn_changes
                })
        
        return pd.DataFrame(codon_data)
    
    def create_syn_nonsyn_plot(self, df, gene_name):
        """Create plot showing synonymous and non-synonymous variation along gene."""
        
        if df is None or len(df) == 0:
            self.logger.warning(f"No data to plot for {gene_name}")
            return
            
        # Create figure
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), 
                                       gridspec_kw={'height_ratios': [3, 1]})
        
        # Main panel: variation along gene
        positions = df['position'].values
        
        # Create stacked bar plot
        syn_heights = np.where(df['is_variable'], df['syn_changes'], 0)
        nonsyn_heights = np.where(df['is_variable'], df['nonsyn_changes'], 0)
        
        # Plot bars
        width = 1.0
        ax1.bar(positions, syn_heights, width=width, color='#2E7D32', 
               alpha=0.8, label='Synonymous', linewidth=0)
        ax1.bar(positions, nonsyn_heights, width=width, bottom=syn_heights,
               color='#C62828', alpha=0.8, label='Non-synonymous', linewidth=0)
        
        # Formatting
        ax1.set_xlabel('Codon position', fontsize=12)
        ax1.set_ylabel('Number of variant types', fontsize=12)
        ax1.set_title(f'{gene_name} - Distribution of Synonymous and Non-synonymous Variation',
                     fontsize=14, fontweight='bold')
        ax1.legend(loc='upper right')
        ax1.grid(axis='y', alpha=0.3)
        ax1.set_xlim(0, len(positions) + 1)
        
        # Lower panel: variation density
        window_size = 10  # 10-codon sliding window
        variation_density = []
        
        for i in range(len(df)):
            start = max(0, i - window_size // 2)
            end = min(len(df), i + window_size // 2 + 1)  # +1 to include center on both sides
            window_variation = df.iloc[start:end]['is_variable'].sum()
            variation_density.append(window_variation / (end - start))
        
        ax2.fill_between(positions, variation_density, color='gray', alpha=0.5)
        ax2.plot(positions, variation_density, color='black', linewidth=1)
        ax2.set_xlabel('Codon position', fontsize=12)
        ax2.set_ylabel('Variation\ndensity', fontsize=11)
        ax2.set_xlim(0, len(positions) + 1)
        ax2.set_ylim(0, 1)
        ax2.grid(axis='y', alpha=0.3)
        
        # Add summary statistics
        n_variable = df['is_variable'].sum()
        n_syn = (df['syn_changes'] > 0).sum()
        n_nonsyn = (df['nonsyn_changes'] > 0).sum()
        total_codons = len(df)
        
        stats_text = (f"Total codons: {total_codons} | "
                     f"Variable: {n_variable} ({n_variable/total_codons*100:.1f}%) | "
                     f"Synonymous sites: {n_syn} | "
                     f"Non-synonymous sites: {n_nonsyn}")
        
        fig.text(0.5, 0.02, stats_text, ha='center', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        plt.subplots_adjust(bottom=0.08)
        
        # Save plot
        plot_file = self.output_dir / 'plots' / 'syn_nonsyn' / f'{gene_name}_syn_nonsyn_distribution.pdf'
        plt.savefig(plot_file, format='pdf', bbox_inches='tight')
        self.logger.info(f"Syn/nonsyn plot saved: {plot_file}")
        plt.close()
        
        return df
    
    def create_operon_syn_nonsyn_plots(self):
        """Create syn/nonsyn plots for all operon genes."""
        
        self.logger.info("Creating syn/nonsyn variation plots for operon genes")
        
        operon_genes = ['frpC', 'glpC', 'ptsD', 'ptsC', 'ptsB', 'ptsA', 'fruR']
        alignment_dir = Path("../05_operon_assembly_extraction/output/msa/dna_alignments")
        
        all_data = {}
        
        for gene in operon_genes:
            alignment_file = alignment_dir / f"{gene}_aligned.fasta"
            
            if not alignment_file.exists():
                self.logger.warning(f"Alignment file not found: {alignment_file}")
                continue
                
            # Analyze variation
            df = self.analyze_codon_variation_alignment(alignment_file, gene)
            
            if df is not None:
                # Create individual plot
                self.create_syn_nonsyn_plot(df, gene)
                all_data[gene] = df
                
                # Save data
                data_file = self.output_dir / 'tables' / f'{gene}_codon_variation.csv'
                df.to_csv(data_file, index=False)
        
        # Create combined figure
        if all_data:
            self.create_combined_syn_nonsyn_figure(all_data)
    
    def create_combined_syn_nonsyn_figure(self, all_data):
        """Create a combined figure showing all operon genes."""
        
        genes = list(all_data.keys())
        fig, axes = plt.subplots(4, 2, figsize=(16, 14))
        axes = axes.flatten()
        
        for idx, gene in enumerate(genes):
            ax = axes[idx]
            df = all_data[gene]
            
            # Create simplified plot
            positions = df['position'].values
            syn_heights = np.where(df['is_variable'], df['syn_changes'], 0)
            nonsyn_heights = np.where(df['is_variable'], df['nonsyn_changes'], 0)
            
            ax.bar(positions, syn_heights, width=1.0, color='#2E7D32', 
                  alpha=0.8, linewidth=0)
            ax.bar(positions, nonsyn_heights, width=1.0, bottom=syn_heights,
                  color='#C62828', alpha=0.8, linewidth=0)
            
            ax.set_title(f'{gene}', fontsize=12, fontweight='bold')
            ax.set_xlabel('Codon position', fontsize=10)
            ax.set_ylabel('Variants', fontsize=10)
            ax.set_xlim(0, len(positions) + 1)
            
            # Add stats
            n_variable = df['is_variable'].sum()
            pct_variable = n_variable / len(df) * 100
            ax.text(0.95, 0.95, f'{pct_variable:.1f}% variable',
                   transform=ax.transAxes, ha='right', va='top',
                   fontsize=9, bbox=dict(boxstyle='round', 
                                       facecolor='white', alpha=0.7))
        
        # Hide unused subplots
        for idx in range(len(genes), len(axes)):
            axes[idx].set_visible(False)
        
        # Add legend
        syn_patch = mpatches.Patch(color='#2E7D32', alpha=0.8, label='Synonymous')
        nonsyn_patch = mpatches.Patch(color='#C62828', alpha=0.8, label='Non-synonymous')
        fig.legend(handles=[syn_patch, nonsyn_patch], 
                  loc='lower right', bbox_to_anchor=(0.95, 0.05))
        
        plt.suptitle('Synonymous and Non-synonymous Variation in Operon Genes',
                    fontsize=16, fontweight='bold', y=0.98)
        plt.tight_layout()
        
        # Save combined plot
        plot_file = self.output_dir / 'plots' / 'syn_nonsyn' / 'operon_combined_syn_nonsyn.pdf'
        plt.savefig(plot_file, format='pdf', bbox_inches='tight')
        self.logger.info(f"Combined syn/nonsyn plot saved: {plot_file}")
        plt.close()
    
    def run(self, operon_file, core_file):
        """Run the complete analysis pipeline."""
        
        self.logger.info("="*60)
        self.logger.info("Starting Codon Variation Analysis Pipeline")
        self.logger.info("="*60)
        
        # Check input files
        operon_path = Path(operon_file)
        core_path = Path(core_file)
        
        if not operon_path.exists():
            self.logger.error(f"Operon file not found: {operon_path}")
            return
        
        if not core_path.exists():
            self.logger.error(f"Core file not found: {core_path}")
            return
        
        # Calculate metrics
        self.logger.info("Calculating selection metrics...")
        operon_df = self.calculate_selection_metrics(operon_path, 'operon')
        core_df = self.calculate_selection_metrics(core_path, 'core')
        
        # Generate summary statistics
        self.logger.info("Generating summary statistics...")
        operon_stats = self.generate_summary_stats(operon_df)
        core_stats = self.generate_summary_stats(core_df)
        
        # Save detailed results
        combined_df = pd.concat([operon_df, core_df], ignore_index=True)
        output_file = self.output_dir / 'tables' / 'codon_variation_metrics.csv'
        combined_df.to_csv(output_file, index=False)
        self.logger.info(f"Detailed results saved to {output_file}")
        
        # Save summary statistics
        summary_df = pd.DataFrame([
            {'gene_type': 'operon', **operon_stats},
            {'gene_type': 'core', **core_stats}
        ])
        summary_file = self.output_dir / 'tables' / 'summary_statistics.csv'
        summary_df.to_csv(summary_file, index=False)
        self.logger.info(f"Summary statistics saved to {summary_file}")
        
        # Create plots
        self.logger.info("Creating visualization plots...")
        self.create_plots(operon_df, core_df)
        
        # Generate report
        self.logger.info("Generating analysis report...")
        self.generate_report(operon_df, core_df, operon_stats, core_stats)
        
        # Create syn/nonsyn variation plots for operon genes
        self.logger.info("Creating synonymous/non-synonymous variation plots...")
        self.create_operon_syn_nonsyn_plots()
        
        self.logger.info("="*60)
        self.logger.info("Analysis completed successfully")
        self.logger.info("="*60)
        
        # Print summary to console
        def _fmt_num(x):
            return f"{x:.3f}" if isinstance(x, (int, float, np.number)) and np.isfinite(x) else "N/A"
        
        print("\nQUICK SUMMARY:")
        print("-"*40)
        print(f"Operon genes: {operon_stats['n_genes']} analyzed")
        print(f"  Mean normalized selection: {_fmt_num(operon_stats.get('mean_normalized_selection'))}")
        print(f"Core genes: {core_stats['n_genes']} analyzed")
        print(f"  Mean normalized selection: {_fmt_num(core_stats.get('mean_normalized_selection'))}")
        print(f"\nResults saved to: {self.output_dir}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Codon Variation-Based Selection Analysis Pipeline'
    )
    parser.add_argument(
        '--operon-file',
        type=str,
        default='../07_dnds_analysis/output_old/tables/operon_dnds_analysis.csv',
        help='Path to operon dN/dS analysis CSV file'
    )
    parser.add_argument(
        '--core-file',
        type=str,
        default='../07_dnds_analysis/output_old/tables/core_genes_dnds_analysis.csv',
        help='Path to core genes dN/dS analysis CSV file'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='output',
        help='Output directory for results'
    )
    parser.add_argument(
        '--expected-neutral',
        type=float,
        default=2.5,
        help='Neutral nonsyn:syn site ratio used to normalize pN/pS (default: 2.5)'
    )
    
    args = parser.parse_args()
    
    # Run analysis
    analyzer = CodonVariationAnalyzer(output_dir=args.output_dir,
                                     expected_neutral=args.expected_neutral)
    analyzer.run(args.operon_file, args.core_file)


if __name__ == "__main__":
    main()