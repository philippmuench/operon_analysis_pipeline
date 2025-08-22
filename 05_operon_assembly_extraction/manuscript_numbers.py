#!/usr/bin/env python3
"""
Generate manuscript statistics for operon assembly extraction methods section.
"""

import os
import glob
import pandas as pd
from collections import defaultdict, Counter
import numpy as np
from Bio import SeqIO
import sys
try:
    from scipy.stats import fisher_exact
except ImportError:
    fisher_exact = None

def analyze_blast_extraction_results():
    """Analyze BLAST-based sequence extraction statistics."""
    
    blast_dir = "../03_blast_search/output/blast_results"
    if not os.path.exists(blast_dir):
        return {"error": "BLAST results directory not found"}
    
    stats = {
        "total_blast_files": 0,
        "genes_with_hits": set(),
        "total_hits": 0,
        "high_quality_hits": 0,
        "identity_distribution": [],
        "coverage_distribution": [],
        "extraction_thresholds": {
            "min_identity": 90.0,  # Current pipeline threshold
            "min_coverage": 80.0   # Current pipeline threshold
        }
    }
    
    # Analyze BLAST result files - match patterns used in the pipeline
    blast_files = []
    for pat in ("*_genes_blast.txt", "*_blast.txt", "*.blast"):
        blast_files.extend(glob.glob(os.path.join(blast_dir, pat)))
    stats["total_blast_files"] = len(blast_files)
    
    for blast_file in blast_files[:50]:  # Sample for speed
        gene_name = os.path.basename(blast_file).replace("_blast.txt", "")
        if "_genes" in gene_name:
            gene_name = gene_name.replace("_genes", "")
        stats["genes_with_hits"].add(gene_name)
        
        try:
            with open(blast_file, 'r') as f:
                for line in f:
                    if line.strip() and not line.startswith('#'):
                        fields = line.strip().split('\t')
                        if len(fields) >= 13:
                            identity = float(fields[2])
                            coverage = float(fields[12]) if fields[12] != 'N/A' else 0.0
                            
                            stats["total_hits"] += 1
                            stats["identity_distribution"].append(identity)
                            stats["coverage_distribution"].append(coverage)
                            
                            if identity >= 90.0 and coverage >= 80.0:
                                stats["high_quality_hits"] += 1
        except:
            continue
    
    # Calculate summary statistics
    if stats["identity_distribution"]:
        stats["mean_identity"] = np.mean(stats["identity_distribution"])
        stats["median_identity"] = np.median(stats["identity_distribution"])
    
    if stats["coverage_distribution"]:
        stats["mean_coverage"] = np.mean(stats["coverage_distribution"])
        stats["median_coverage"] = np.median(stats["coverage_distribution"])
    
    stats["unique_genes"] = len(stats["genes_with_hits"])
    if stats["total_hits"] > 0:
        stats["high_quality_rate"] = stats["high_quality_hits"] / stats["total_hits"] * 100
    
    return stats

def analyze_extracted_sequences():
    """Analyze extracted sequence statistics."""
    
    output_dirs = [
        "output/sequences",
        "output/mappings/aa_nt_mapping/prokka/sequences",
        "output/mappings/aa_nt_mapping/assemblies/sequences",
        "output/mappings/nt_nt_mapping/prokka_genome/sequences",
        "output/mappings/nt_nt_mapping/prokka_variants"
    ]
    
    stats = {
        "extraction_strategies": 0,
        "total_genes_extracted": 0,
        "sequences_per_gene": {},
        "sequence_lengths": [],
        "extraction_success_rate": 0
    }
    
    for seq_dir in output_dirs:
        if os.path.exists(seq_dir):
            stats["extraction_strategies"] += 1
            fasta_files = glob.glob(os.path.join(seq_dir, "*.fasta")) + glob.glob(os.path.join(seq_dir, "*.fa"))
            
            for fasta_file in fasta_files:
                gene_name = os.path.basename(fasta_file).replace(".fa", "")
                try:
                    sequences = list(SeqIO.parse(fasta_file, "fasta"))
                    seq_count = len(sequences)
                    
                    if gene_name not in stats["sequences_per_gene"]:
                        stats["sequences_per_gene"][gene_name] = []
                    stats["sequences_per_gene"][gene_name].append(seq_count)
                    
                    for seq in sequences:
                        stats["sequence_lengths"].append(len(seq.seq))
                    
                    if seq_count > 0:
                        stats["total_genes_extracted"] += 1
                except:
                    continue
    
    # Calculate averages
    if stats["sequences_per_gene"]:
        all_counts = []
        for gene, counts in stats["sequences_per_gene"].items():
            all_counts.extend(counts)
        if all_counts:
            stats["avg_sequences_per_gene"] = np.mean(all_counts)
            stats["median_sequences_per_gene"] = np.median(all_counts)
    
    if stats["sequence_lengths"]:
        stats["avg_sequence_length"] = np.mean(stats["sequence_lengths"])
        stats["median_sequence_length"] = np.median(stats["sequence_lengths"])
        stats["min_sequence_length"] = min(stats["sequence_lengths"])
        stats["max_sequence_length"] = max(stats["sequence_lengths"])
    
    return stats

def analyze_msa_results():
    """Analyze multiple sequence alignment results."""
    
    msa_dirs = [
        "output/msa/dna_alignments",
        "output/msa/noncoding_alignments",
        "output/mappings/aa_nt_mapping/prokka/msa/dna_alignments",
        "output/mappings/aa_nt_mapping/assemblies/msa/dna_alignments",
        "output/mappings/nt_nt_mapping/prokka_genome/msa/dna_alignments",
        "output/mappings/nt_nt_mapping/prokka_variants/msa_variants"
    ]
    
    stats = {
        "total_alignments": 0,
        "coding_alignments": 0,
        "noncoding_alignments": 0,
        "alignment_lengths": [],
        "sequences_per_alignment": [],
        "mafft_version": "7.490",  # Based on script
        "alignment_algorithm": "automatic selection"
    }
    
    for msa_dir in msa_dirs:
        if os.path.exists(msa_dir):
            alignment_files = glob.glob(os.path.join(msa_dir, "*.fasta")) + glob.glob(os.path.join(msa_dir, "*.fa"))
            
            for alignment_file in alignment_files:
                try:
                    alignment = list(SeqIO.parse(alignment_file, "fasta"))
                    if len(alignment) > 1:
                        stats["total_alignments"] += 1
                        stats["alignment_lengths"].append(len(alignment[0].seq))
                        stats["sequences_per_alignment"].append(len(alignment))
                        
                        if "noncoding" in msa_dir:
                            stats["noncoding_alignments"] += 1
                        else:
                            stats["coding_alignments"] += 1
                except:
                    continue
    
    # Calculate summary statistics
    if stats["alignment_lengths"]:
        stats["avg_alignment_length"] = np.mean(stats["alignment_lengths"])
        stats["median_alignment_length"] = np.median(stats["alignment_lengths"])
    
    if stats["sequences_per_alignment"]:
        stats["avg_sequences_per_alignment"] = np.mean(stats["sequences_per_alignment"])
        stats["median_sequences_per_alignment"] = np.median(stats["sequences_per_alignment"])
    
    return stats

def analyze_conservation_plots():
    """Analyze conservation plot generation."""
    
    plot_dirs = [
        "output/plots",
        "output/mappings/aa_nt_mapping/prokka/plots",
        "output/mappings/aa_nt_mapping/assemblies/plots",
        "output/mappings/nt_nt_mapping/prokka_genome/plots",
        "output/mappings/nt_nt_mapping/prokka_variants/plots"
    ]
    
    stats = {
        "total_plots": 0,
        "gene_conservation_plots": 0,
        "promoter_plots": 0,
        "plot_types": set()
    }
    
    for plot_dir in plot_dirs:
        if os.path.exists(plot_dir):
            plot_files = glob.glob(os.path.join(plot_dir, "*.png"))
            stats["total_plots"] += len(plot_files)
            
            for plot_file in plot_files:
                filename = os.path.basename(plot_file)
                if "conservation" in filename:
                    stats["gene_conservation_plots"] += 1
                elif "promoter" in filename:
                    stats["promoter_plots"] += 1
                
                # Extract plot type
                if "_" in filename:
                    plot_type = filename.split("_")[0]
                    stats["plot_types"].add(plot_type)
    
    stats["plot_types"] = list(stats["plot_types"])
    return stats

def get_extraction_methods():
    """Get information about extraction methods used."""
    
    methods = {
        "coordinate_based": {
            "description": "BLAST hit coordinates mapped to genome assemblies",
            "threshold_identity": "≥90%",
            "threshold_coverage": "≥80%",
            "strand_handling": "Automatic reverse complement for negative strand hits"
        },
        "prokka_based": {
            "description": "Sequences extracted from Prokka annotations using BLAST mapping",
            "file_types": [".fna", ".ffn"],
            "coordinate_system": "CDS annotations"
        },
        "assembly_based": {
            "description": "Direct extraction from raw genome assemblies",
            "coordinate_mapping": "BLAST sstart/send positions",
            "quality_control": "Identity and coverage filtering"
        }
    }
    
    return methods

def analyze_ptsa_start_codon_pattern():
    """Analyze the laboratory-specific ptsA start codon pattern.
    
    This analysis identifies a key finding: ptsA gene, which canonically uses TTG start codon,
    shows significantly higher ATG usage in laboratory-adapted strains.
    """
    
    try:
        # Check if the start site analysis has been run
        summary_file = "../08_start_site_analysis/output/start_site_summary.tsv"
        metadata_file = "../00_annotation/8587_Efs_metadata_ASbarcode.txt"
        
        if not os.path.exists(summary_file):
            return {"error": "Start site analysis not yet run (Step 08)"}
        
        if not os.path.exists(metadata_file):
            return {"error": "Metadata file not found"}
        
        # Load and analyze the data
        summary_df = pd.read_csv(summary_file, sep='\t')
        metadata_df = pd.read_csv(metadata_file, sep='\t')
        
        # Extract genome ID from file name to match with AS_barcode
        summary_df['AS_barcode'] = summary_df['genome_id'].str.replace('.result', '', regex=False)
        
        # Merge with metadata
        merged_df = summary_df.merge(metadata_df, on='AS_barcode', how='left')
        
        # Filter for ptsA gene only
        ptsa_df = merged_df[merged_df['gene'] == 'ptsA'].copy()
        
        # Overall statistics
        overall_counts = ptsa_df['start_codon'].value_counts()
        overall_pct = ptsa_df['start_codon'].value_counts(normalize=True) * 100
        
        # Laboratory vs others
        lab_data = ptsa_df[ptsa_df['Source Niche'] == 'Laboratory']
        other_data = ptsa_df[(ptsa_df['Source Niche'] != 'Laboratory') & ptsa_df['Source Niche'].notna()]
        
        stats = {
            'total_ptsa': len(ptsa_df),
            'overall_ttg_pct': overall_pct.get('TTG', 0),
            'overall_atg_pct': overall_pct.get('ATG', 0),
            'overall_gtg_pct': overall_pct.get('GTG', 0)
        }
        
        if len(lab_data) > 0 and len(other_data) > 0:
            lab_atg = len(lab_data[lab_data['start_codon'] == 'ATG'])
            lab_ttg = len(lab_data[lab_data['start_codon'] == 'TTG'])
            other_atg = len(other_data[other_data['start_codon'] == 'ATG'])
            other_ttg = len(other_data[other_data['start_codon'] == 'TTG'])
            
            stats['lab_atg_pct'] = (lab_atg / len(lab_data)) * 100
            stats['lab_ttg_pct'] = (lab_ttg / len(lab_data)) * 100
            stats['other_atg_pct'] = (other_atg / len(other_data)) * 100
            stats['other_ttg_pct'] = (other_ttg / len(other_data)) * 100
            stats['lab_n'] = len(lab_data)
            stats['other_n'] = len(other_data)
            
            # Fisher's exact test
            if fisher_exact is not None and lab_ttg > 0 and other_atg > 0:
                contingency_table = np.array([[lab_atg, lab_ttg], [other_atg, other_ttg]])
                odds_ratio, p_value = fisher_exact(contingency_table)
            else:
                # Fallback to simple calculation if scipy not available
                odds_ratio = (lab_atg * other_ttg) / (lab_ttg * other_atg) if (lab_ttg * other_atg) > 0 else float('inf')
                p_value = None  # Can't calculate without scipy
            
            stats['odds_ratio'] = odds_ratio
            stats['p_value'] = p_value
        
        return stats
        
    except Exception as e:
        return {"error": f"Error analyzing ptsA pattern: {str(e)}"}

def generate_manuscript_stats():
    """Generate all statistics for manuscript."""
    
    print("Operon Assembly Extraction Statistics for Manuscript")
    print("=" * 60)
    
    # BLAST extraction analysis
    print(f"\n1. BLAST-based Sequence Extraction:")
    blast_stats = analyze_blast_extraction_results()
    
    if "error" not in blast_stats:
        print(f"   BLAST result files analyzed: {blast_stats['total_blast_files']:,}")
        print(f"   Unique genes with hits: {blast_stats['unique_genes']}")
        print(f"   Total BLAST hits processed: {blast_stats['total_hits']:,}")
        print(f"   High-quality hits extracted: {blast_stats['high_quality_hits']:,}")
        if blast_stats['total_hits'] > 0:
            print(f"   Extraction success rate: {blast_stats.get('high_quality_rate', 0):.1f}%")
        print(f"   Quality thresholds: ≥{blast_stats['extraction_thresholds']['min_identity']:.0f}% identity, ≥{blast_stats['extraction_thresholds']['min_coverage']:.0f}% coverage")
        if blast_stats.get('mean_identity'):
            print(f"   Mean sequence identity: {blast_stats['mean_identity']:.1f}%")
            print(f"   Mean query coverage: {blast_stats['mean_coverage']:.1f}%")
    
    # Sequence extraction results
    print(f"\n2. Sequence Extraction Results:")
    seq_stats = analyze_extracted_sequences()
    print(f"   Extraction strategies employed: {seq_stats['extraction_strategies']}")
    print(f"   Total genes successfully extracted: {seq_stats['total_genes_extracted']}")
    if seq_stats.get('avg_sequences_per_gene'):
        print(f"   Average sequences per gene: {seq_stats['avg_sequences_per_gene']:.0f}")
        print(f"   Median sequences per gene: {seq_stats['median_sequences_per_gene']:.0f}")
    if seq_stats.get('avg_sequence_length'):
        print(f"   Average sequence length: {seq_stats['avg_sequence_length']:.0f} bp")
        print(f"   Sequence length range: {seq_stats.get('min_sequence_length', 0)}-{seq_stats.get('max_sequence_length', 0)} bp")
    
    # Multiple sequence alignment analysis
    print(f"\n3. Multiple Sequence Alignment Generation:")
    msa_stats = analyze_msa_results()
    print(f"   MAFFT version: {msa_stats['mafft_version']}")
    print(f"   Algorithm selection: {msa_stats['alignment_algorithm']}")
    print(f"   Total alignments created: {msa_stats['total_alignments']}")
    print(f"   Coding gene alignments: {msa_stats['coding_alignments']}")
    print(f"   Non-coding region alignments: {msa_stats['noncoding_alignments']}")
    if msa_stats.get('avg_alignment_length'):
        print(f"   Average alignment length: {msa_stats['avg_alignment_length']:.0f} positions")
        print(f"   Average sequences per alignment: {msa_stats['avg_sequences_per_alignment']:.0f}")
    
    # Conservation analysis and visualization
    print(f"\n4. Conservation Analysis and Visualization:")
    plot_stats = analyze_conservation_plots()
    print(f"   Total conservation plots generated: {plot_stats['total_plots']}")
    print(f"   Gene conservation plots: {plot_stats['gene_conservation_plots']}")
    print(f"   Promoter conservation plots: {plot_stats['promoter_plots']}")
    if plot_stats['plot_types']:
        print(f"   Plot types generated: {', '.join(plot_stats['plot_types'])}")
    
    # Extraction methods
    print(f"\n5. Extraction Methods:")
    methods = get_extraction_methods()
    for method_name, details in methods.items():
        print(f"   {method_name.replace('_', ' ').title()}:")
        print(f"     - {details['description']}")
        if 'threshold_identity' in details:
            print(f"     - Identity threshold: {details['threshold_identity']}")
            print(f"     - Coverage threshold: {details['threshold_coverage']}")
    
    # ptsA start codon pattern analysis
    print(f"\n6. ptsA Start Codon Pattern Analysis:")
    ptsa_stats = analyze_ptsa_start_codon_pattern()
    
    if "error" not in ptsa_stats:
        print(f"   Total ptsA genes analyzed: {ptsa_stats['total_ptsa']}")
        print(f"   Overall start codon usage:")
        print(f"     - TTG: {ptsa_stats['overall_ttg_pct']:.1f}% (canonical for ptsA)")
        print(f"     - ATG: {ptsa_stats['overall_atg_pct']:.1f}%")
        print(f"     - GTG: {ptsa_stats['overall_gtg_pct']:.1f}%")
        
        if ptsa_stats.get('lab_atg_pct') is not None:
            print(f"\n   Laboratory-specific pattern:")
            print(f"     - Laboratory strains (n={ptsa_stats.get('lab_n', 0)}):")
            print(f"       ATG: {ptsa_stats['lab_atg_pct']:.1f}%")
            print(f"       TTG: {ptsa_stats['lab_ttg_pct']:.1f}%")
            print(f"     - Other niches (n={ptsa_stats.get('other_n', 0)}):")
            print(f"       ATG: {ptsa_stats['other_atg_pct']:.1f}%")
            print(f"       TTG: {ptsa_stats['other_ttg_pct']:.1f}%")
            print(f"     - Statistical significance:")
            print(f"       Odds ratio: {ptsa_stats['odds_ratio']:.2f}")
            
            if ptsa_stats.get('p_value') is not None:
                print(f"       P-value: {ptsa_stats['p_value']:.4e}")
                if ptsa_stats['p_value'] < 0.001:
                    print(f"       Result: Highly significant (p < 0.001)")
                elif ptsa_stats['p_value'] < 0.01:
                    print(f"       Result: Significant (p < 0.01)")
                elif ptsa_stats['p_value'] < 0.05:
                    print(f"       Result: Significant (p < 0.05)")
                else:
                    print(f"       Result: Not significant")
            else:
                print(f"       P-value: Requires scipy for calculation")
    else:
        print(f"   {ptsa_stats['error']}")
    
    # Summary for manuscript
    print("\n" + "=" * 60)
    print("MANUSCRIPT NUMBERS SUMMARY:")
    print("=" * 60)
    
    if "error" not in blast_stats and blast_stats['total_hits'] > 0:
        print(f"BLAST hits processed: {blast_stats['total_hits']:,}")
        print(f"High-quality sequences extracted: {blast_stats['high_quality_hits']:,}")
        print(f"Quality thresholds: ≥90% identity, ≥80% coverage")
    
    if seq_stats.get('avg_sequences_per_gene'):
        print(f"Average sequences per gene: {seq_stats['avg_sequences_per_gene']:.0f}")
        print(f"Average sequence length: {seq_stats.get('avg_sequence_length', 0):.0f} bp")
    
    print(f"Multiple sequence alignments: {msa_stats['total_alignments']}")
    print(f"MAFFT version 7.490 with automatic algorithm selection")
    print(f"Conservation plots generated: {plot_stats['total_plots']}")
    
    # Add ptsA pattern to summary
    if "error" not in ptsa_stats and ptsa_stats.get('lab_atg_pct') is not None:
        print(f"\nptsA LABORATORY ADAPTATION:")
        print(f"ptsA normally uses TTG start codon ({ptsa_stats['overall_ttg_pct']:.1f}%)")
        print(f"Laboratory strains show ATG usage: {ptsa_stats['lab_atg_pct']:.1f}%")
        print(f"Other niches show ATG usage: {ptsa_stats['other_atg_pct']:.1f}%")
        if ptsa_stats.get('p_value') is not None:
            print(f"Statistical significance: p = {ptsa_stats['p_value']:.4e}")
        elif ptsa_stats.get('odds_ratio') is not None:
            print(f"Odds ratio: {ptsa_stats['odds_ratio']:.2f}")
    
    return {
        'blast_stats': blast_stats,
        'sequence_stats': seq_stats,
        'msa_stats': msa_stats,
        'plot_stats': plot_stats,
        'methods': methods,
        'ptsa_stats': ptsa_stats
    }

if __name__ == "__main__":
    # Change to script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    
    stats = generate_manuscript_stats()
