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
from math import log
from itertools import combinations
from functools import lru_cache
import warnings
warnings.filterwarnings('ignore')

# Static cached function for syn_fraction_by_pos
@lru_cache(maxsize=4096)
def _syn_fraction_by_pos_cached_static(codon, table_id=11):
    """Static cached version of syn_fraction_by_pos."""
    NUCS = "ATCG"
    try:
        aa0 = str(Seq(codon).translate(table=table_id))
    except:
        return None
        
    if aa0 == "*" or any(n not in NUCS for n in codon):  # skip stops/ambiguous
        return None
    
    fracs = []
    for pos in range(3):
        syn = 0
        tot = 0
        for n in NUCS:
            if n == codon[pos]:
                continue
            mutated = codon[:pos] + n + codon[pos+1:]
            if any(x not in NUCS for x in mutated):
                continue
            try:
                aa = str(Seq(mutated).translate(table=table_id))
            except:
                continue
            if aa == "*":
                continue  # exclude paths through stop
            tot += 1
            if aa == aa0:
                syn += 1
        fracs.append(syn / tot if tot else 0.0)
    return tuple(fracs)

class DNdSAnalyzer:
    """Main class for dN/dS analysis pipeline."""
    
    def __init__(self, operon_dir, core_dir, output_dir, threads=1, include_simple=False):
        self.operon_dir = Path(operon_dir) if operon_dir else None
        self.core_dir = Path(core_dir) if core_dir else None
        self.output_dir = Path(output_dir)
        self.threads = threads
        self.include_simple = include_simple
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
        
    def classify_substitution_simple(self, codon1, codon2):
        """Simple classification of codon substitution (original method for comparison)."""
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
    
    def syn_fraction_by_pos_cached(self, codon, table_id=11):
        """Cached version of syn_fraction_by_pos."""
        return _syn_fraction_by_pos_cached_static(codon, table_id)
    
    def syn_fraction_by_pos(self, codon, table_id=11):
        """For each position (0..2) return fraction of single-nuc changes that are synonymous."""
        NUCS = "ATCG"
        try:
            aa0 = str(Seq(codon).translate(table=table_id))
        except:
            return None
            
        if aa0 == "*" or any(n not in NUCS for n in codon):  # skip stops/ambiguous
            return None
        
        fracs = []
        for pos in range(3):
            syn = 0
            tot = 0
            for n in NUCS:
                if n == codon[pos]:
                    continue
                mutated = codon[:pos] + n + codon[pos+1:]
                if any(x not in NUCS for x in mutated):
                    continue
                try:
                    aa = str(Seq(mutated).translate(table=table_id))
                except:
                    continue
                if aa == "*":
                    continue  # exclude paths through stop
                tot += 1
                if aa == aa0:
                    syn += 1
            fracs.append(syn / tot if tot else 0.0)
        return fracs
    
    def pairwise_dn_ds_ng86(self, codons_A, codons_B, table_id=11):
        """NG86-style dN/dS calculation between two sequences' codon strings.
        
        Implements proper site normalization and JC correction.
        """
        NUCS = "ATCG"
        S_sites = N_sites = 0.0
        syn_diffs = nonsyn_diffs = 0.0

        for c1, c2 in zip(codons_A, codons_B):
            if "-" in c1 or "-" in c2 or len(c1) != 3 or len(c2) != 3:
                continue
            if any(n not in NUCS for n in c1 + c2):
                continue
            
            # Let the cached function handle stop codons and ambiguous bases
            f1 = self.syn_fraction_by_pos_cached(c1, table_id)
            f2 = self.syn_fraction_by_pos_cached(c2, table_id)
            if f1 is None or f2 is None:
                continue  # covers stops/ambigs
            
            # Average site opportunities across the two codons
            for p in range(3):
                fS = 0.5 * (f1[p] + f2[p])
                S_sites += fS
                N_sites += (1.0 - fS)

            # Distribute observed diffs across S/N using those fractions
            for p in range(3):
                if c1[p] != c2[p]:
                    fS = 0.5 * (f1[p] + f2[p])
                    syn_diffs += fS
                    nonsyn_diffs += (1.0 - fS)

        # Convert to proportions
        pS = syn_diffs / S_sites if S_sites > 0 else float("nan")
        pN = nonsyn_diffs / N_sites if N_sites > 0 else float("nan")

        # JC correction (Nei–Gojobori)
        def jc_correction(p):
            if p is None or p != p or p >= 0.75:
                return float("inf") if p >= 0.75 else float("nan")
            return -0.75 * log(1 - (4.0/3.0) * p)

        dS = jc_correction(pS)
        dN = jc_correction(pN)
        omega = dN / dS if (dS > 0 and dS != float("inf") and not (dS != dS)) else float("nan")
        
        return dN, dS, omega, S_sites, N_sites, syn_diffs, nonsyn_diffs
    
    def reference_based_dn_ds(self, alignment, gene_name=None, table_id=11):
        """Calculate dN/dS using reference-based approach for large alignments.
        
        Uses actual reference sequence from step 02 when available.
        Much faster for alignments with many sequences.
        """
        from collections import Counter
        
        # Try to load reference sequence
        ref_file = Path("../02_reference_operon_extraction/output/operon_genes_nt.fasta")
        reference_seq = None
        
        if gene_name and ref_file.exists():
            try:
                for record in SeqIO.parse(ref_file, "fasta"):
                    # Parse header: >operon_CDS_1|frpC|fructoselysine-6-phosphate_deglycase
                    parts = record.id.split('|')
                    if len(parts) >= 2 and parts[1] == gene_name:
                        reference_seq = str(record.seq).upper()
                        self.logger.info(f"Using reference sequence for {gene_name}")
                        break
            except Exception as e:
                self.logger.warning(f"Could not load reference: {e}")
        
        # Build consensus per codon position to preserve frame
        seq_len = alignment.get_alignment_length()
        ref_codons = []
        
        if reference_seq:
            # Use provided reference
            ref_codons = [reference_seq[i:i+3] for i in range(0, len(reference_seq)-2, 3)]
        else:
            # Build consensus codon-wise from alignment positions
            for i in range(0, seq_len, 3):
                if i+2 >= seq_len:
                    break
                # Get columns for this codon
                cols = [alignment[:, i+k] for k in range(3)]
                # Majority vote per position
                codon = ''
                for col in cols:
                    # Count bases excluding gaps
                    bases = [b for b in col if b in 'ACGT']
                    if bases:
                        # Use most common non-gap base
                        codon += Counter(bases).most_common(1)[0][0]
                    else:
                        # If no valid bases, mark as N
                        codon += 'N'
                ref_codons.append(codon)
            self.logger.info(f"Using consensus sequence as reference for {gene_name or 'unknown'}")
        
        # Precompute per-position synonymous fractions for reference
        fS_by_pos = []
        for rc in ref_codons:
            if len(rc) != 3 or '-' in rc or 'N' in rc:
                fS_by_pos.append(None)
                continue
            fr = self.syn_fraction_by_pos_cached(rc, table_id)
            fS_by_pos.append(fr)
        
        # Calculate sites and differences
        total_S_sites = 0.0
        total_N_sites = 0.0
        syn_diffs_total = 0.0
        nonsyn_diffs_total = 0.0
        n_seqs = 0
        
        for record in alignment:
            s = str(record.seq).upper()
            # Extract codons from the same positions as reference
            seq_codons = [s[i:i+3] for i in range(0, seq_len, 3) if i+2 < seq_len]
            n_seqs += 1
            
            for rc, sc, fr in zip(ref_codons, seq_codons, fS_by_pos):
                if fr is None:
                    continue
                # Only skip if BOTH reference and sequence have invalid codons
                if len(rc) != 3 or len(sc) != 3:
                    continue
                # Skip if either has N or if both have gaps
                if 'N' in rc or 'N' in sc:
                    continue
                if '-' in rc and '-' in sc:
                    continue
                # Allow comparison if only one has gaps (treat as difference)
                if '-' in rc or '-' in sc:
                    # Can't properly calculate syn/nonsyn for gap comparisons
                    continue
                    
                # Add site opportunities (only once per sequence)
                for p in range(3):
                    total_S_sites += fr[p]
                    total_N_sites += (1.0 - fr[p])
                    
                    # If position differs, add fractional diffs
                    if rc[p] != sc[p]:
                        syn_diffs_total += fr[p]
                        nonsyn_diffs_total += (1.0 - fr[p])
        
        # Calculate proportions
        pS = syn_diffs_total / total_S_sites if total_S_sites > 0 else float('nan')
        pN = nonsyn_diffs_total / total_N_sites if total_N_sites > 0 else float('nan')
        
        # Apply Jukes-Cantor correction
        def jc_correction(p):
            if p <= 0:
                return 0.0
            if p >= 0.749999:
                return float('inf')
            return -0.75 * log(1 - 4.0*p/3.0)
        
        dS = jc_correction(pS)
        dN = jc_correction(pN)
        
        omega = dN / dS if (np.isfinite(dS) and dS > 0) else (float('inf') if dN > 0 else float('nan'))
        
        return dN, dS, omega, total_S_sites, total_N_sites, syn_diffs_total, nonsyn_diffs_total
            
    def analyze_alignment(self, alignment_file):
        """Analyze a single alignment for dN/dS and SNP types."""
        try:
            alignment = AlignIO.read(alignment_file, 'fasta')
            n_seqs = len(alignment)
            seq_length = alignment.get_alignment_length()
            
            if n_seqs < 2:
                return None
                
            # Extract gene name from filename
            gene_name = os.path.basename(alignment_file).replace('_aligned.fasta', '').replace('.fasta', '')
            
            # Initialize counters
            stats = {
                'file': os.path.basename(alignment_file),
                'gene': gene_name,
                'n_sequences': n_seqs,
                'alignment_length': seq_length,
                'variable_sites': 0,
                'parsimony_informative': 0,
                'singleton_sites': 0,
                'gap_percentage': 0
            }
            
            # Simple method (optional)
            if self.include_simple:
                stats['synonymous'] = 0
                stats['non_synonymous'] = 0
                stats['total_codons'] = 0
                
                # Analyze codons
                for pos in range(0, seq_length - 2, 3):
                    codons = [str(seq[pos:pos+3].seq).upper() for seq in alignment]
                    codons = [c for c in codons if '-' not in c and len(c) == 3]
                    
                    if len(codons) < 2:
                        continue
                        
                    stats['total_codons'] += len(codons)
                    
                    # Pairwise comparisons (simple method)
                    for i in range(len(codons)-1):
                        for j in range(i+1, len(codons)):
                            change_type, _ = self.classify_substitution_simple(codons[i], codons[j])
                            if change_type == 'synonymous':
                                stats['synonymous'] += 1
                            elif change_type == 'non-synonymous':
                                stats['non_synonymous'] += 1
                                
                # Calculate simple dN/dS (no site normalization)
                if stats['synonymous'] > 0:
                    stats['simple_dn_ds'] = stats['non_synonymous'] / stats['synonymous']
                else:
                    stats['simple_dn_ds'] = np.nan if stats['non_synonymous'] == 0 else np.inf
            
            # Method 2: NG86 with proper site normalization and JC correction
            # Build codon strings for each sequence
            seq_codons = []
            total_valid_codons = 0
            for record in alignment:
                s = str(record.seq).upper()
                # Skip sequences not divisible by 3 or with too many gaps
                if len(s) % 3 == 0:
                    codons = [s[i:i+3] for i in range(0, len(s), 3) if i+3 <= len(s)]
                    valid_codons = [c for c in codons if '-' not in c and 'N' not in c and len(c) == 3]
                    total_valid_codons += len(valid_codons)
                    seq_codons.append(codons)
            
            # Add total valid codons count
            stats['ng86_total_codons'] = total_valid_codons
            
            
            # For large alignments (>100 sequences), use reference-based approach
            if n_seqs > 100:
                self.logger.info(f"Using reference-based approach for {gene_name} ({n_seqs} sequences)")
                dN, dS, omega, S_sites, N_sites, syn_diffs, nonsyn_diffs = self.reference_based_dn_ds(
                    alignment, gene_name
                )
                
                # Calculate total valid codons for reference-based approach
                # S_sites + N_sites already represents the sum across all sequences
                # Divide by 3 to convert from nucleotide sites to codons
                total_codons_ref = (S_sites + N_sites) / 3.0
                
                # Add NG86 statistics from reference-based calculation
                stats['ng86_total_codons'] = int(round(total_codons_ref))
                stats['ng86_dn_ds'] = omega
                stats['ng86_mean_dn_ds'] = omega
                stats['ng86_mean_dN'] = dN
                stats['ng86_mean_dS'] = dS
                stats['ng86_S_sites'] = S_sites
                stats['ng86_N_sites'] = N_sites
                stats['ng86_syn_diffs'] = syn_diffs
                stats['ng86_nonsyn_diffs'] = nonsyn_diffs
                stats['ng86_n_pairs'] = n_seqs  # Number of comparisons to reference
                stats['ng86_method'] = 'reference_based'
            else:
                # For small alignments, use pairwise comparisons
                self.logger.info(f"Using pairwise approach for {gene_name} ({n_seqs} sequences)")
                
                # Aggregate rates instead of collecting individual omegas
                sum_dN = 0.0
                sum_dS = 0.0
                sum_S_sites = 0.0
                sum_N_sites = 0.0
                pairs_used = 0
                
                if len(seq_codons) >= 2:
                    for i, j in combinations(range(len(seq_codons)), 2):
                        dN, dS, omega, S_sites, N_sites, syn_diffs, nonsyn_diffs = self.pairwise_dn_ds_ng86(
                            seq_codons[i], seq_codons[j]
                        )
                        
                        # Aggregate the rates
                        if np.isfinite(dN):
                            sum_dN += dN
                        if np.isfinite(dS):
                            sum_dS += dS
                        if np.isfinite(S_sites):
                            sum_S_sites += S_sites
                        if np.isfinite(N_sites):
                            sum_N_sites += N_sites
                        if np.isfinite(dN) or np.isfinite(dS):
                            pairs_used += 1
                
                # Calculate aggregated statistics
                if sum_dS > 0:
                    stats['ng86_dn_ds'] = sum_dN / sum_dS
                elif sum_dN > 0:
                    stats['ng86_dn_ds'] = float('inf')
                else:
                    stats['ng86_dn_ds'] = np.nan
                    
                stats['ng86_mean_dn_ds'] = stats['ng86_dn_ds']  # For compatibility
                stats['ng86_mean_dN'] = sum_dN / max(1, pairs_used)
                stats['ng86_mean_dS'] = sum_dS / max(1, pairs_used)
                stats['ng86_S_sites'] = sum_S_sites / max(1, pairs_used)
                stats['ng86_N_sites'] = sum_N_sites / max(1, pairs_used)
                stats['ng86_n_pairs'] = pairs_used
                stats['ng86_method'] = 'pairwise'
            
            # Keep dn_ds_ratio for backward compatibility, but use NG86 value
            stats['dn_ds_ratio'] = stats['ng86_dn_ds']
            
            # --- Codon-level SNP-type counts (site-based, simple & readable) ---
            codon_variable_sites = 0
            codon_syn_sites = 0
            codon_nonsyn_sites = 0

            for pos in range(0, seq_length - 2, 3):
                observed = []
                for rec in alignment:
                    c = str(rec.seq[pos:pos+3]).upper()
                    if len(c) != 3 or '-' in c or 'N' in c:
                        continue
                    try:
                        aa = str(Seq(c).translate(table=11))
                    except Exception:
                        continue
                    if aa == '*':  # skip alleles translating to stop
                        continue
                    observed.append((c, aa))

                # need at least 2 valid codons to be variable
                if len(observed) < 2:
                    continue

                codon_set = {c for c, _ in observed}
                if len(codon_set) < 2:
                    continue  # not variable at codon level

                aa_set = {aa for _, aa in observed}
                codon_variable_sites += 1
                if len(aa_set) == 1:
                    codon_syn_sites += 1
                else:
                    codon_nonsyn_sites += 1

            stats['codon_variable_sites'] = codon_variable_sites
            stats['codon_syn_sites'] = codon_syn_sites
            stats['codon_nonsyn_sites'] = codon_nonsyn_sites
            stats['codon_syn_fraction'] = (codon_syn_sites / codon_variable_sites
                                           if codon_variable_sites > 0 else np.nan)
            # --- end codon-level SNP counts ---
                
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
        
    def create_method_comparison_plot(self, df):
        """Create plot comparing simple vs NG86 methods."""
        # Only create comparison if simple method was computed
        if 'simple_dn_ds' not in df.columns or df['simple_dn_ds'].isna().all():
            self.logger.info("Skipping method comparison plot (simple method not computed)")
            return
        
        if 'ng86_dn_ds' not in df.columns:
            return
        
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        fig.suptitle('Comparison of dN/dS Methods', fontsize=14, fontweight='bold')
        
        # Filter finite values
        simple_vals = df['simple_dn_ds'].replace([np.inf, -np.inf], np.nan).dropna()
        ng86_vals = df['ng86_dn_ds'].replace([np.inf, -np.inf], np.nan).dropna()
        
        # Scatter plot comparing methods
        ax1 = axes[0]
        merged = df[['gene', 'simple_dn_ds', 'ng86_dn_ds']].dropna()
        merged = merged[np.isfinite(merged['simple_dn_ds']) & np.isfinite(merged['ng86_dn_ds'])]
        if not merged.empty:
            ax1.scatter(merged['simple_dn_ds'], merged['ng86_dn_ds'], alpha=0.6)
            ax1.plot([0, 2], [0, 2], 'r--', alpha=0.5, label='y=x')
            ax1.set_xlabel('Simple dN/dS (no normalization)')
            ax1.set_ylabel('NG86 dN/dS (site-normalized + JC)')
            ax1.set_title('Method Comparison')
            ax1.legend()
            ax1.grid(True, alpha=0.3)
            
            # Add gene labels for outliers
            for _, row in merged.iterrows():
                if abs(row['simple_dn_ds'] - row['ng86_dn_ds']) > 0.5:
                    ax1.annotate(row['gene'], (row['simple_dn_ds'], row['ng86_dn_ds']), 
                               fontsize=8, alpha=0.7)
        
        # Distribution comparison
        ax2 = axes[1]
        ax2.hist(simple_vals, bins=20, alpha=0.5, label='Simple method', color='blue')
        ax2.hist(ng86_vals, bins=20, alpha=0.5, label='NG86 method', color='red')
        ax2.set_xlabel('dN/dS ratio')
        ax2.set_ylabel('Frequency')
        ax2.set_title('Distribution Comparison')
        ax2.legend()
        ax2.set_xlim(0, 2)
        
        # Boxplot comparison
        ax3 = axes[2]
        plot_data = []
        labels = []
        if len(simple_vals) > 0:
            plot_data.append(simple_vals)
            labels.append(f'Simple\n(mean={simple_vals.mean():.3f})')
        if len(ng86_vals) > 0:
            plot_data.append(ng86_vals)
            labels.append(f'NG86\n(mean={ng86_vals.mean():.3f})')
        
        if plot_data:
            bp = ax3.boxplot(plot_data, labels=labels)
            ax3.set_ylabel('dN/dS ratio')
            ax3.set_title('Method Comparison')
            ax3.axhline(y=1, color='r', linestyle='--', alpha=0.5, label='Neutral evolution')
            ax3.legend()
        
        plt.tight_layout()
        plot_file = self.output_dir / 'plots' / 'method_comparison.png'
        plt.savefig(plot_file, dpi=150, bbox_inches='tight')
        plt.close()
        self.logger.info(f"Method comparison plot saved to {plot_file}")
    
    def create_visualizations(self, operon_df, core_df):
        """Create comprehensive visualizations."""
        # First create method comparison if both methods are available
        if not operon_df.empty:
            self.create_method_comparison_plot(operon_df)
        
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        fig.suptitle('dN/dS Analysis Results (NG86 Method)', fontsize=16, fontweight='bold')
        
        # 1. dN/dS distribution comparison (using NG86 method)
        ax = axes[0, 0]
        # Use NG86 values if available, otherwise fall back to simple
        ratio_col = 'ng86_dn_ds' if 'ng86_dn_ds' in operon_df.columns else 'dn_ds_ratio'
        if not operon_df.empty and ratio_col in operon_df.columns:
            operon_ratios = operon_df[ratio_col].dropna()
            operon_ratios = operon_ratios[np.isfinite(operon_ratios)]
            ax.hist(operon_ratios, bins=30, alpha=0.7, 
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
                
                # Report NG86 method results (primary)
                if 'ng86_dn_ds' in operon_df.columns:
                    f.write(f"\n=== NG86 METHOD (Site-normalized with JC correction) ===\n")
                    ratios = operon_df['ng86_dn_ds'].dropna()
                    ratios_finite = ratios[np.isfinite(ratios)]
                    
                    f.write(f"\ndN/dS Statistics (NG86):\n")
                    if len(ratios_finite) > 0:
                        f.write(f"  Mean: {ratios_finite.mean():.4f}\n")
                        f.write(f"  Median: {ratios_finite.median():.4f}\n")
                        f.write(f"  Std Dev: {ratios_finite.std():.4f}\n")
                        f.write(f"  Min: {ratios_finite.min():.4f}\n")
                        f.write(f"  Max: {ratios_finite.max():.4f}\n")
                        f.write(f"  25th percentile: {ratios_finite.quantile(0.25):.4f}\n")
                        f.write(f"  75th percentile: {ratios_finite.quantile(0.75):.4f}\n")
                        
                        f.write(f"\nSelection Categories (NG86):\n")
                        f.write(f"  Strong purifying (ω < 0.1): {sum(ratios_finite < 0.1)} genes\n")
                        f.write(f"  Purifying (0.1 ≤ ω < 0.5): {sum((ratios_finite >= 0.1) & (ratios_finite < 0.5))} genes\n")
                        f.write(f"  Weak purifying (0.5 ≤ ω < 1): {sum((ratios_finite >= 0.5) & (ratios_finite < 1))} genes\n")
                        f.write(f"  Neutral/positive (ω ≥ 1): {sum(ratios_finite >= 1)} genes\n")
                    
                    if 'ng86_mean_dN' in operon_df.columns:
                        f.write(f"\nRate Statistics (NG86):\n")
                        f.write(f"  Mean dN: {operon_df['ng86_mean_dN'].mean():.6f}\n")
                        f.write(f"  Mean dS: {operon_df['ng86_mean_dS'].mean():.6f}\n")
                        f.write(f"  Mean S-sites: {operon_df['ng86_S_sites'].mean():.1f}\n")
                        f.write(f"  Mean N-sites: {operon_df['ng86_N_sites'].mean():.1f}\n")
                    
                    if 'ng86_syn_diffs' in operon_df.columns:
                        f.write(f"\nSubstitution Counts (NG86 method):\n")
                        f.write(f"  Mean synonymous differences: {operon_df['ng86_syn_diffs'].mean():.2f}\n")
                        f.write(f"  Mean non-synonymous differences: {operon_df['ng86_nonsyn_diffs'].mean():.2f}\n")
                        f.write(f"  Total valid codons analyzed: {operon_df['ng86_total_codons'].sum()}\n")
                
                # Report site-based SNP counts
                if 'codon_variable_sites' in operon_df.columns:
                        f.write(f"\nSite-based SNP Counts (simple, readable):\n")
                        f.write(f"  Total variable codon sites: {operon_df['codon_variable_sites'].sum()}\n")
                        f.write(f"  Synonymous sites: {operon_df['codon_syn_sites'].sum()}\n")
                        f.write(f"  Non-synonymous sites: {operon_df['codon_nonsyn_sites'].sum()}\n")
                        f.write(f"  Mean synonymous fraction: {operon_df['codon_syn_fraction'].mean():.3f}\n")
                
                # Report simple method for comparison
                if 'simple_dn_ds' in operon_df.columns:
                    f.write(f"\n=== SIMPLE METHOD (Count-based, no normalization) ===\n")
                    simple_ratios = operon_df['simple_dn_ds'].dropna()
                    simple_finite = simple_ratios[np.isfinite(simple_ratios)]
                    
                    if len(simple_finite) > 0:
                        f.write(f"\ndN/dS Statistics (Simple):\n")
                        f.write(f"  Mean: {simple_finite.mean():.4f}\n")
                        f.write(f"  Median: {simple_finite.median():.4f}\n")
                        f.write(f"  Note: Simple method overestimates ω due to lack of site normalization\n")
                
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
        
        # Log analysis mode
        if self.operon_dir and self.core_dir:
            self.logger.info("Mode: Analyzing both operon and core genes")
        elif self.operon_dir:
            self.logger.info("Mode: Analyzing operon genes only")
        elif self.core_dir:
            self.logger.info("Mode: Analyzing core genes only")
        
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
    parser.add_argument(
        '--include-simple',
        action='store_true',
        help='Also compute the naive count ratio (slower for large datasets)'
    )
    parser.add_argument(
        '--mode',
        type=str,
        choices=['both', 'operon', 'core'],
        default='both',
        help='Which gene sets to analyze: both (default), operon, or core'
    )
    
    args = parser.parse_args()
    
    # Handle mode-based directory requirements
    if args.mode == 'operon':
        if not args.operon_dir:
            # Use default operon directory
            args.operon_dir = '../05_operon_assembly_extraction/output/msa/dna_alignments'
        args.core_dir = None  # Disable core analysis
    elif args.mode == 'core':
        if not args.core_dir:
            # Use default core directory
            args.core_dir = '../04_core_gene_analysis/output/core_gene_alignments'
        args.operon_dir = None  # Disable operon analysis
    else:  # mode == 'both'
        # Use defaults if not provided
        if not args.operon_dir:
            args.operon_dir = '../05_operon_assembly_extraction/output/msa/dna_alignments'
        if not args.core_dir:
            args.core_dir = '../04_core_gene_analysis/output/core_gene_alignments'
    
    # Validate that at least one input directory is provided
    if not args.operon_dir and not args.core_dir:
        parser.error("No input directories available for the selected mode")
    
    # Run the pipeline
    analyzer = DNdSAnalyzer(
        operon_dir=args.operon_dir,
        core_dir=args.core_dir,
        output_dir=args.output_dir,
        threads=args.threads,
        include_simple=args.include_simple
    )
    
    analyzer.run()


if __name__ == "__main__":
    main()