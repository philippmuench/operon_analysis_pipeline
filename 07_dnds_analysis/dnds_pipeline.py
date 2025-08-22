#!/usr/bin/env python3
"""
Consolidated dN/dS Analysis Pipeline for Operon and Core Genes.
Combines all analysis steps into a single comprehensive pipeline.
"""

import os
import sys
import glob
import argparse
import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from collections import defaultdict, Counter
from multiprocessing import Pool, cpu_count
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

class DNdSAnalyzer:
    """Main class for dN/dS analysis pipeline."""
    
    def __init__(self, operon_dir, core_dir, output_dir, threads=1):
        self.operon_dir = Path(operon_dir) if operon_dir else None
        self.core_dir = Path(core_dir) if core_dir else None
        self.output_dir = Path(output_dir)
        self.threads = threads
        self.genetic_code = CodonTable.unambiguous_dna_by_id[11]
        
        # Create output directory structure
        self.output_dir.mkdir(parents=True, exist_ok=True)
        (self.output_dir / 'plots').mkdir(exist_ok=True)
        (self.output_dir / 'tables').mkdir(exist_ok=True)
        
        # Set up logging
        self.setup_logging()
        
    def setup_logging(self):
        """Set up logging configuration."""
        log_file = self.output_dir / 'dnds_pipeline.log'
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
        
    def classify_substitution(self, codon1, codon2):
        """Classify a codon substitution as synonymous or non-synonymous."""
        if codon1 == codon2:
            return 'identical', 0
            
        if len(codon1) != 3 or len(codon2) != 3:
            return 'invalid', 0
            
        # Handle ambiguous nucleotides
        if any(n not in 'ATCG' for n in codon1 + codon2):
            return 'ambiguous', 0
            
        try:
            aa1 = str(Seq(codon1).translate(table=11))
            aa2 = str(Seq(codon2).translate(table=11))
            
            # Count nucleotide differences
            n_diffs = sum(1 for i in range(3) if codon1[i] != codon2[i])
            
            if aa1 == aa2:
                return 'synonymous', n_diffs
            else:
                return 'non-synonymous', n_diffs
        except:
            return 'invalid', 0
            
    def analyze_alignment(self, alignment_file):
        """Analyze a single alignment for dN/dS and SNP types."""
        try:
            alignment = AlignIO.read(alignment_file, 'fasta')
            n_seqs = len(alignment)
            seq_length = alignment.get_alignment_length()
            
            if n_seqs < 2:
                return None
                
            # Initialize counters
            stats = {
                'file': os.path.basename(alignment_file),
                'n_sequences': n_seqs,
                'alignment_length': seq_length,
                'synonymous': 0,
                'non_synonymous': 0,
                'total_codons': 0,
                'variable_sites': 0,
                'parsimony_informative': 0,
                'singleton_sites': 0,
                'gap_percentage': 0
            }
            
            # Analyze codons
            for pos in range(0, seq_length - 2, 3):
                codons = [str(seq[pos:pos+3].seq).upper() for seq in alignment]
                codons = [c for c in codons if '-' not in c and len(c) == 3]
                
                if len(codons) < 2:
                    continue
                    
                stats['total_codons'] += len(codons)
                
                # Pairwise comparisons
                for i in range(len(codons)-1):
                    for j in range(i+1, len(codons)):
                        change_type, _ = self.classify_substitution(codons[i], codons[j])
                        if change_type == 'synonymous':
                            stats['synonymous'] += 1
                        elif change_type == 'non-synonymous':
                            stats['non_synonymous'] += 1
                            
            # Calculate dN/dS
            if stats['synonymous'] > 0:
                stats['dn_ds_ratio'] = stats['non_synonymous'] / stats['synonymous']
            else:
                stats['dn_ds_ratio'] = np.nan if stats['non_synonymous'] == 0 else np.inf
                
            # Calculate site statistics
            for pos in range(seq_length):
                column = alignment[:, pos]
                chars = [c for c in column if c not in '-N']
                
                if len(chars) < 2:
                    continue
                    
                char_counts = Counter(chars)
                if len(char_counts) > 1:
                    stats['variable_sites'] += 1
                    
                    # Check if parsimony informative (at least 2 chars occur >= 2 times)
                    if sum(1 for count in char_counts.values() if count >= 2) >= 2:
                        stats['parsimony_informative'] += 1
                    elif len(char_counts) == 2 and min(char_counts.values()) == 1:
                        stats['singleton_sites'] += 1
                        
            # Calculate gap percentage
            total_chars = n_seqs * seq_length
            gap_count = sum(str(seq.seq).count('-') for seq in alignment)
            stats['gap_percentage'] = (gap_count / total_chars) * 100
            
            return stats
            
        except Exception as e:
            self.logger.error(f"Error analyzing {alignment_file}: {e}")
            return None
            
    def process_alignments(self, alignment_dir, pattern="*.fasta"):
        """Process all alignments in a directory."""
        alignment_files = list(alignment_dir.glob(pattern))
        if not alignment_files:
            alignment_files = list(alignment_dir.glob("*_aligned.fasta"))
            
        self.logger.info(f"Found {len(alignment_files)} alignment files in {alignment_dir}")
        
        if self.threads > 1:
            with Pool(processes=self.threads) as pool:
                results = pool.map(self.analyze_alignment, alignment_files)
        else:
            results = [self.analyze_alignment(f) for f in alignment_files]
            
        # Filter out None results
        results = [r for r in results if r is not None]
        
        return pd.DataFrame(results) if results else pd.DataFrame()
        
    def create_visualizations(self, operon_df, core_df):
        """Create comprehensive visualizations."""
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle('dN/dS Analysis Results', fontsize=16, fontweight='bold')
        
        # 1. dN/dS distribution comparison
        ax = axes[0, 0]
        if not operon_df.empty and 'dn_ds_ratio' in operon_df.columns:
            operon_ratios = operon_df['dn_ds_ratio'].dropna()
            ax.hist(operon_ratios[operon_ratios != np.inf], bins=30, alpha=0.7, 
                   label=f'Operon genes (n={len(operon_ratios)})', color='blue')
        
        if not core_df.empty and 'dn_ds_ratio' in core_df.columns:
            core_ratios = core_df['dn_ds_ratio'].dropna()
            ax.hist(core_ratios[core_ratios != np.inf], bins=30, alpha=0.7,
                   label=f'Core genes (n={len(core_ratios)})', color='red')
        
        ax.axvline(x=1, color='black', linestyle='--', alpha=0.5, label='Neutral (dN/dS=1)')
        ax.set_xlabel('dN/dS ratio')
        ax.set_ylabel('Frequency')
        ax.set_title('Distribution of dN/dS Ratios')
        ax.legend()
        ax.set_xlim(0, 3)
        
        # 2. Substitution types
        ax = axes[0, 1]
        categories = ['Operon', 'Core']
        syn_means = []
        nonsyn_means = []
        
        if not operon_df.empty:
            syn_means.append(operon_df['synonymous'].mean() if 'synonymous' in operon_df else 0)
            nonsyn_means.append(operon_df['non_synonymous'].mean() if 'non_synonymous' in operon_df else 0)
        else:
            syn_means.append(0)
            nonsyn_means.append(0)
            
        if not core_df.empty:
            syn_means.append(core_df['synonymous'].mean() if 'synonymous' in core_df else 0)
            nonsyn_means.append(core_df['non_synonymous'].mean() if 'non_synonymous' in core_df else 0)
        else:
            syn_means.append(0)
            nonsyn_means.append(0)
            
        x = np.arange(len(categories))
        width = 0.35
        
        ax.bar(x - width/2, syn_means, width, label='Synonymous', color='green', alpha=0.7)
        ax.bar(x + width/2, nonsyn_means, width, label='Non-synonymous', color='orange', alpha=0.7)
        ax.set_xlabel('Gene Category')
        ax.set_ylabel('Mean Substitutions')
        ax.set_title('Average Substitution Types')
        ax.set_xticks(x)
        ax.set_xticklabels(categories)
        ax.legend()
        
        # 3. Selection pressure categories
        ax = axes[0, 2]
        selection_categories = {
            'Strong purifying (<0.1)': lambda x: x < 0.1,
            'Purifying (0.1-0.5)': lambda x: (x >= 0.1) & (x < 0.5),
            'Weak purifying (0.5-1)': lambda x: (x >= 0.5) & (x < 1),
            'Neutral/positive (≥1)': lambda x: x >= 1
        }
        
        operon_counts = []
        core_counts = []
        
        for category, condition in selection_categories.items():
            if not operon_df.empty and 'dn_ds_ratio' in operon_df.columns:
                operon_ratios = operon_df['dn_ds_ratio'].dropna()
                operon_ratios = operon_ratios[operon_ratios != np.inf]
                operon_counts.append(sum(condition(operon_ratios)))
            else:
                operon_counts.append(0)
                
            if not core_df.empty and 'dn_ds_ratio' in core_df.columns:
                core_ratios = core_df['dn_ds_ratio'].dropna()
                core_ratios = core_ratios[core_ratios != np.inf]
                core_counts.append(sum(condition(core_ratios)))
            else:
                core_counts.append(0)
        
        x = np.arange(len(selection_categories))
        width = 0.35
        
        ax.bar(x - width/2, operon_counts, width, label='Operon', color='blue', alpha=0.7)
        ax.bar(x + width/2, core_counts, width, label='Core', color='red', alpha=0.7)
        ax.set_xlabel('Selection Category')
        ax.set_ylabel('Number of Genes')
        ax.set_title('Selection Pressure Categories')
        ax.set_xticks(x)
        ax.set_xticklabels([cat.split('(')[0].strip() for cat in selection_categories.keys()], 
                          rotation=45, ha='right')
        ax.legend()
        
        # 4. Variable sites analysis
        ax = axes[1, 0]
        if not operon_df.empty and 'variable_sites' in operon_df.columns:
            ax.scatter(operon_df['alignment_length'], operon_df['variable_sites'], 
                      alpha=0.6, label='Operon', color='blue')
        if not core_df.empty and 'variable_sites' in core_df.columns:
            ax.scatter(core_df['alignment_length'], core_df['variable_sites'], 
                      alpha=0.6, label='Core', color='red')
        ax.set_xlabel('Alignment Length (bp)')
        ax.set_ylabel('Variable Sites')
        ax.set_title('Variable Sites vs Alignment Length')
        ax.legend()
        
        # 5. Gap percentage distribution
        ax = axes[1, 1]
        if not operon_df.empty and 'gap_percentage' in operon_df.columns:
            ax.hist(operon_df['gap_percentage'], bins=20, alpha=0.7, 
                   label='Operon', color='blue')
        if not core_df.empty and 'gap_percentage' in core_df.columns:
            ax.hist(core_df['gap_percentage'], bins=20, alpha=0.7,
                   label='Core', color='red')
        ax.set_xlabel('Gap Percentage (%)')
        ax.set_ylabel('Frequency')
        ax.set_title('Distribution of Gap Percentages')
        ax.legend()
        
        # 6. Statistical summary box
        ax = axes[1, 2]
        ax.axis('off')
        
        summary_text = "Statistical Summary\n" + "="*30 + "\n\n"
        
        if not operon_df.empty and 'dn_ds_ratio' in operon_df.columns:
            operon_ratios = operon_df['dn_ds_ratio'].dropna()
            operon_ratios = operon_ratios[operon_ratios != np.inf]
            summary_text += f"Operon Genes:\n"
            summary_text += f"  Mean dN/dS: {operon_ratios.mean():.3f}\n"
            summary_text += f"  Median dN/dS: {operon_ratios.median():.3f}\n"
            summary_text += f"  Std dN/dS: {operon_ratios.std():.3f}\n"
            summary_text += f"  Under selection: {sum(operon_ratios < 1)}/{len(operon_ratios)}\n\n"
        
        if not core_df.empty and 'dn_ds_ratio' in core_df.columns:
            core_ratios = core_df['dn_ds_ratio'].dropna()
            core_ratios = core_ratios[core_ratios != np.inf]
            summary_text += f"Core Genes:\n"
            summary_text += f"  Mean dN/dS: {core_ratios.mean():.3f}\n"
            summary_text += f"  Median dN/dS: {core_ratios.median():.3f}\n"
            summary_text += f"  Std dN/dS: {core_ratios.std():.3f}\n"
            summary_text += f"  Under selection: {sum(core_ratios < 1)}/{len(core_ratios)}\n\n"
        
        # Statistical test if both datasets available
        if (not operon_df.empty and 'dn_ds_ratio' in operon_df.columns and 
            not core_df.empty and 'dn_ds_ratio' in core_df.columns):
            operon_ratios = operon_df['dn_ds_ratio'].dropna()
            operon_ratios = operon_ratios[operon_ratios != np.inf]
            core_ratios = core_df['dn_ds_ratio'].dropna()
            core_ratios = core_ratios[core_ratios != np.inf]
            
            if len(operon_ratios) > 0 and len(core_ratios) > 0:
                statistic, pvalue = stats.mannwhitneyu(operon_ratios, core_ratios, alternative='two-sided')
                summary_text += f"Mann-Whitney U test:\n"
                summary_text += f"  p-value: {pvalue:.2e}\n"
                summary_text += f"  Significant: {'Yes' if pvalue < 0.05 else 'No'}"
        
        ax.text(0.1, 0.9, summary_text, transform=ax.transAxes, 
               fontsize=10, verticalalignment='top', fontfamily='monospace')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'plots' / 'dnds_analysis_summary.png', dpi=300, bbox_inches='tight')
        plt.savefig(self.output_dir / 'plots' / 'dnds_analysis_summary.pdf', dpi=300, bbox_inches='tight')
        plt.close()
        
        self.logger.info(f"Saved visualization to {self.output_dir / 'plots'}")
        
    def generate_summary_report(self, operon_df, core_df):
        """Generate a comprehensive text summary report."""
        report_file = self.output_dir / 'dnds_summary_report.txt'
        
        with open(report_file, 'w') as f:
            f.write("="*80 + "\n")
            f.write("dN/dS ANALYSIS SUMMARY REPORT\n")
            f.write("="*80 + "\n\n")
            
            # Operon genes summary
            f.write("OPERON GENES ANALYSIS\n")
            f.write("-"*40 + "\n")
            
            if not operon_df.empty:
                f.write(f"Total genes analyzed: {len(operon_df)}\n")
                
                if 'dn_ds_ratio' in operon_df.columns:
                    ratios = operon_df['dn_ds_ratio'].dropna()
                    ratios_finite = ratios[ratios != np.inf]
                    
                    f.write(f"\ndN/dS Statistics:\n")
                    f.write(f"  Mean: {ratios_finite.mean():.4f}\n")
                    f.write(f"  Median: {ratios_finite.median():.4f}\n")
                    f.write(f"  Std Dev: {ratios_finite.std():.4f}\n")
                    f.write(f"  Min: {ratios_finite.min():.4f}\n")
                    f.write(f"  Max: {ratios_finite.max():.4f}\n")
                    f.write(f"  25th percentile: {ratios_finite.quantile(0.25):.4f}\n")
                    f.write(f"  75th percentile: {ratios_finite.quantile(0.75):.4f}\n")
                    
                    f.write(f"\nSelection Categories:\n")
                    f.write(f"  Strong purifying (dN/dS < 0.1): {sum(ratios_finite < 0.1)}\n")
                    f.write(f"  Purifying (0.1 ≤ dN/dS < 0.5): {sum((ratios_finite >= 0.1) & (ratios_finite < 0.5))}\n")
                    f.write(f"  Weak purifying (0.5 ≤ dN/dS < 1): {sum((ratios_finite >= 0.5) & (ratios_finite < 1))}\n")
                    f.write(f"  Neutral/positive (dN/dS ≥ 1): {sum(ratios_finite >= 1)}\n")
                
                if 'synonymous' in operon_df.columns and 'non_synonymous' in operon_df.columns:
                    f.write(f"\nSubstitution Statistics:\n")
                    f.write(f"  Total synonymous: {operon_df['synonymous'].sum()}\n")
                    f.write(f"  Total non-synonymous: {operon_df['non_synonymous'].sum()}\n")
                    f.write(f"  Mean synonymous per gene: {operon_df['synonymous'].mean():.2f}\n")
                    f.write(f"  Mean non-synonymous per gene: {operon_df['non_synonymous'].mean():.2f}\n")
                
                if 'variable_sites' in operon_df.columns:
                    f.write(f"\nSequence Variation:\n")
                    f.write(f"  Mean variable sites: {operon_df['variable_sites'].mean():.2f}\n")
                    f.write(f"  Mean parsimony informative sites: {operon_df['parsimony_informative'].mean():.2f}\n")
                    f.write(f"  Mean singleton sites: {operon_df['singleton_sites'].mean():.2f}\n")
                    f.write(f"  Mean gap percentage: {operon_df['gap_percentage'].mean():.2f}%\n")
                
                # List individual genes
                f.write(f"\nIndividual Gene Results:\n")
                for _, row in operon_df.iterrows():
                    if 'dn_ds_ratio' in row:
                        f.write(f"  {row['file']}: dN/dS = {row['dn_ds_ratio']:.4f}\n")
            else:
                f.write("No operon gene alignments found or analyzed.\n")
            
            # Core genes summary
            f.write("\n" + "="*80 + "\n")
            f.write("CORE GENES ANALYSIS\n")
            f.write("-"*40 + "\n")
            
            if not core_df.empty:
                f.write(f"Total genes analyzed: {len(core_df)}\n")
                
                if 'dn_ds_ratio' in core_df.columns:
                    ratios = core_df['dn_ds_ratio'].dropna()
                    ratios_finite = ratios[ratios != np.inf]
                    
                    f.write(f"\ndN/dS Statistics:\n")
                    f.write(f"  Mean: {ratios_finite.mean():.4f}\n")
                    f.write(f"  Median: {ratios_finite.median():.4f}\n")
                    f.write(f"  Std Dev: {ratios_finite.std():.4f}\n")
                    f.write(f"  Min: {ratios_finite.min():.4f}\n")
                    f.write(f"  Max: {ratios_finite.max():.4f}\n")
                    f.write(f"  25th percentile: {ratios_finite.quantile(0.25):.4f}\n")
                    f.write(f"  75th percentile: {ratios_finite.quantile(0.75):.4f}\n")
                    
                    f.write(f"\nSelection Categories:\n")
                    f.write(f"  Strong purifying (dN/dS < 0.1): {sum(ratios_finite < 0.1)}\n")
                    f.write(f"  Purifying (0.1 ≤ dN/dS < 0.5): {sum((ratios_finite >= 0.1) & (ratios_finite < 0.5))}\n")
                    f.write(f"  Weak purifying (0.5 ≤ dN/dS < 1): {sum((ratios_finite >= 0.5) & (ratios_finite < 1))}\n")
                    f.write(f"  Neutral/positive (dN/dS ≥ 1): {sum(ratios_finite >= 1)}\n")
                
                if 'synonymous' in core_df.columns and 'non_synonymous' in core_df.columns:
                    f.write(f"\nSubstitution Statistics:\n")
                    f.write(f"  Total synonymous: {core_df['synonymous'].sum()}\n")
                    f.write(f"  Total non-synonymous: {core_df['non_synonymous'].sum()}\n")
                    f.write(f"  Mean synonymous per gene: {core_df['synonymous'].mean():.2f}\n")
                    f.write(f"  Mean non-synonymous per gene: {core_df['non_synonymous'].mean():.2f}\n")
            else:
                f.write("No core gene alignments found or analyzed.\n")
            
            # Comparative analysis
            f.write("\n" + "="*80 + "\n")
            f.write("COMPARATIVE ANALYSIS\n")
            f.write("-"*40 + "\n")
            
            if not operon_df.empty and not core_df.empty:
                if 'dn_ds_ratio' in operon_df.columns and 'dn_ds_ratio' in core_df.columns:
                    operon_ratios = operon_df['dn_ds_ratio'].dropna()
                    operon_ratios = operon_ratios[operon_ratios != np.inf]
                    core_ratios = core_df['dn_ds_ratio'].dropna()
                    core_ratios = core_ratios[core_ratios != np.inf]
                    
                    if len(operon_ratios) > 0 and len(core_ratios) > 0:
                        # Mann-Whitney U test
                        statistic, pvalue = stats.mannwhitneyu(operon_ratios, core_ratios, alternative='two-sided')
                        f.write(f"\nStatistical Test (Mann-Whitney U):\n")
                        f.write(f"  Test statistic: {statistic:.2f}\n")
                        f.write(f"  P-value: {pvalue:.2e}\n")
                        f.write(f"  Significant difference (α=0.05): {'Yes' if pvalue < 0.05 else 'No'}\n")
                        
                        # Effect size (rank-biserial correlation)
                        n1, n2 = len(operon_ratios), len(core_ratios)
                        r = 1 - (2*statistic)/(n1*n2)
                        f.write(f"  Effect size (rank-biserial): {r:.3f}\n")
                        
                        # Interpretation
                        f.write(f"\nInterpretation:\n")
                        mean_diff = operon_ratios.mean() - core_ratios.mean()
                        if pvalue < 0.05:
                            if mean_diff > 0:
                                f.write("  Operon genes show significantly higher dN/dS ratios than core genes,\n")
                                f.write("  suggesting relaxed purifying selection or potential positive selection.\n")
                            else:
                                f.write("  Operon genes show significantly lower dN/dS ratios than core genes,\n")
                                f.write("  suggesting stronger purifying selection.\n")
                        else:
                            f.write("  No significant difference in selection pressure between operon and core genes.\n")
            else:
                f.write("Comparative analysis not possible - missing data for operon or core genes.\n")
            
            f.write("\n" + "="*80 + "\n")
            f.write(f"Report generated: {pd.Timestamp.now()}\n")
            f.write("="*80 + "\n")
        
        self.logger.info(f"Summary report saved to {report_file}")
        
    def run(self):
        """Run the complete dN/dS analysis pipeline."""
        self.logger.info("="*60)
        self.logger.info("Starting dN/dS Analysis Pipeline")
        self.logger.info("="*60)
        
        # Process operon alignments
        operon_df = pd.DataFrame()
        if self.operon_dir and self.operon_dir.exists():
            self.logger.info(f"\nProcessing operon alignments from {self.operon_dir}")
            operon_df = self.process_alignments(self.operon_dir)
            if not operon_df.empty:
                operon_df.to_csv(self.output_dir / 'tables' / 'operon_dnds_analysis.csv', index=False)
                self.logger.info(f"Analyzed {len(operon_df)} operon genes")
        else:
            self.logger.warning(f"Operon directory not found: {self.operon_dir}")
        
        # Process core gene alignments
        core_df = pd.DataFrame()
        if self.core_dir and self.core_dir.exists():
            self.logger.info(f"\nProcessing core gene alignments from {self.core_dir}")
            core_df = self.process_alignments(self.core_dir)
            if not core_df.empty:
                core_df.to_csv(self.output_dir / 'tables' / 'core_genes_dnds_analysis.csv', index=False)
                self.logger.info(f"Analyzed {len(core_df)} core genes")
        else:
            self.logger.warning(f"Core gene directory not found: {self.core_dir}")
        
        # Generate visualizations
        if not operon_df.empty or not core_df.empty:
            self.logger.info("\nGenerating visualizations...")
            self.create_visualizations(operon_df, core_df)
            
            self.logger.info("Generating summary report...")
            self.generate_summary_report(operon_df, core_df)
            
            # Create quick summary statistics file
            summary_stats = {}
            
            if not operon_df.empty and 'dn_ds_ratio' in operon_df.columns:
                ratios = operon_df['dn_ds_ratio'].dropna()
                ratios = ratios[ratios != np.inf]
                summary_stats['operon_mean_dnds'] = ratios.mean()
                summary_stats['operon_median_dnds'] = ratios.median()
                summary_stats['operon_n_genes'] = len(operon_df)
                
            if not core_df.empty and 'dn_ds_ratio' in core_df.columns:
                ratios = core_df['dn_ds_ratio'].dropna()
                ratios = ratios[ratios != np.inf]
                summary_stats['core_mean_dnds'] = ratios.mean()
                summary_stats['core_median_dnds'] = ratios.median()
                summary_stats['core_n_genes'] = len(core_df)
            
            # Save summary statistics
            summary_file = self.output_dir / 'dnds_summary_statistics.txt'
            with open(summary_file, 'w') as f:
                for key, value in summary_stats.items():
                    f.write(f"{key}: {value:.4f}\n" if isinstance(value, float) else f"{key}: {value}\n")
            
            self.logger.info(f"Summary statistics saved to {summary_file}")
        else:
            self.logger.error("No alignment data to analyze")
            
        self.logger.info("\n" + "="*60)
        self.logger.info("Pipeline completed successfully")
        self.logger.info("="*60)
        
        # Print summary to console
        if summary_stats:
            print("\nQuick Summary:")
            print("-"*40)
            for key, value in summary_stats.items():
                print(f"{key}: {value:.4f}" if isinstance(value, float) else f"{key}: {value}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Consolidated dN/dS Analysis Pipeline for Operon and Core Genes'
    )
    parser.add_argument(
        '--operon-dir',
        type=str,
        help='Directory containing operon gene alignments'
    )
    parser.add_argument(
        '--core-dir',
        type=str,
        help='Directory containing core gene alignments'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='output',
        help='Output directory for results (default: output)'
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=1,
        help='Number of threads for parallel processing (default: 1)'
    )
    
    args = parser.parse_args()
    
    # Validate that at least one input directory is provided
    if not args.operon_dir and not args.core_dir:
        parser.error("At least one of --operon-dir or --core-dir must be provided")
    
    # Run the pipeline
    analyzer = DNdSAnalyzer(
        operon_dir=args.operon_dir,
        core_dir=args.core_dir,
        output_dir=args.output_dir,
        threads=args.threads
    )
    
    analyzer.run()


if __name__ == "__main__":
    main()