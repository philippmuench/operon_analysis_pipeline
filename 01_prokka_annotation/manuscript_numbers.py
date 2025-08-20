#!/usr/bin/env python3
"""
Generate manuscript statistics for Prokka annotation methods section.
"""

import os
import glob
import gzip
from collections import defaultdict

def count_input_genomes():
    """Count input genome assemblies."""
    genome_dir = "../../Efs_assemblies"
    if not os.path.exists(genome_dir):
        return 0
    
    genome_files = glob.glob(os.path.join(genome_dir, "*.fasta.gz"))
    return len(genome_files)

def analyze_prokka_outputs():
    """Analyze Prokka annotation outputs with outlier detection."""
    output_dir = "../prokka_output"
    if not os.path.exists(output_dir):
        output_dir = "output/prokka_results"  # Test mode fallback
    
    if not os.path.exists(output_dir):
        return {"error": "Prokka output directory not found"}
    
    # Get all genome directories
    genome_dirs = [d for d in os.listdir(output_dir) 
                   if os.path.isdir(os.path.join(output_dir, d))]
    
    stats = {
        "total_genomes": len(genome_dirs),
        "complete_genomes": 0,
        "total_genes": 0,
        "total_proteins": 0,
        "avg_genes_per_genome": 0,
        "avg_proteins_per_genome": 0,
        "file_counts": defaultdict(int),
        "gene_count_outliers": []
    }
    
    expected_files = ["gff", "faa", "ffn", "fna", "gbk", "log", "txt", "tsv"]
    gene_counts = []
    protein_counts = []
    genome_gene_data = []  # Store (genome_id, gene_count) for outlier detection
    
    # Analyze all genomes (remove sample size limit for comprehensive analysis)
    sample_size = len(genome_dirs)  # Analyze all available genomes
    
    print(f"Starting analysis of {sample_size} genomes...", flush=True)
    progress_interval = max(1, sample_size // 20)  # Show progress 20 times
    
    for idx, genome_dir in enumerate(genome_dirs[:sample_size], 1):
        # Show progress
        if idx % progress_interval == 0 or idx == sample_size:
            print(f"  Progress: {idx}/{sample_size} genomes analyzed ({idx*100//sample_size}%)", flush=True)
        
        genome_path = os.path.join(output_dir, genome_dir)
        
        # Check if genome is complete
        complete = True
        for ext in expected_files:
            file_path = os.path.join(genome_path, f"{genome_dir}.{ext}")
            if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
                stats["file_counts"][ext] += 1
            else:
                complete = False
        
        if complete:
            stats["complete_genomes"] += 1
            
            # Count genes from GFF file
            gff_file = os.path.join(genome_path, f"{genome_dir}.gff")
            if os.path.exists(gff_file):
                try:
                    with open(gff_file, 'r') as f:
                        gene_count = sum(1 for line in f if not line.startswith('#') and 'CDS' in line)
                    gene_counts.append(gene_count)
                    genome_gene_data.append((genome_dir, gene_count))
                    stats["total_genes"] += gene_count
                except:
                    genome_gene_data.append((genome_dir, 0))
            else:
                genome_gene_data.append((genome_dir, 0))
            
            # Count proteins from FAA file
            faa_file = os.path.join(genome_path, f"{genome_dir}.faa")
            if os.path.exists(faa_file):
                try:
                    with open(faa_file, 'r') as f:
                        protein_count = sum(1 for line in f if line.startswith('>'))
                    protein_counts.append(protein_count)
                    stats["total_proteins"] += protein_count
                except:
                    pass
    
    print(f"Analysis complete. Calculating statistics...", flush=True)
    
    # Calculate averages and identify outliers
    if gene_counts:
        import numpy as np
        
        stats["avg_genes_per_genome"] = sum(gene_counts) / len(gene_counts)
        stats["median_genes_per_genome"] = np.median(gene_counts)
        stats["std_genes_per_genome"] = np.std(gene_counts)
        stats["min_genes_per_genome"] = min(gene_counts)
        stats["max_genes_per_genome"] = max(gene_counts)
        
        # Identify outliers (genomes with very few genes)
        # Use Q1 - 1.5*IQR as lower bound for outliers
        q1 = np.percentile(gene_counts, 25)
        q3 = np.percentile(gene_counts, 75)
        iqr = q3 - q1
        lower_bound = q1 - 1.5 * iqr
        
        outliers = []
        very_low_gene_genomes = []
        
        for genome_id, gene_count in genome_gene_data:
            if gene_count < lower_bound:
                outliers.append((genome_id, gene_count))
            if gene_count < 1000:  # Flag genomes with < 1000 genes as potentially problematic
                very_low_gene_genomes.append((genome_id, gene_count))
        
        stats["gene_count_outliers"] = sorted(outliers, key=lambda x: x[1])[:10]  # Top 10 outliers
        stats["very_low_gene_genomes"] = sorted(very_low_gene_genomes, key=lambda x: x[1])[:10]
        stats["num_outliers"] = len(outliers)
        stats["num_very_low_gene"] = len(very_low_gene_genomes)
        stats["lower_bound_threshold"] = lower_bound
        
        # Percentile distribution
        stats["gene_count_percentiles"] = {
            "1%": np.percentile(gene_counts, 1),
            "5%": np.percentile(gene_counts, 5),
            "10%": np.percentile(gene_counts, 10),
            "25%": q1,
            "50%": np.percentile(gene_counts, 50),
            "75%": q3,
            "90%": np.percentile(gene_counts, 90),
            "95%": np.percentile(gene_counts, 95),
            "99%": np.percentile(gene_counts, 99)
        }
    
    if protein_counts:
        stats["avg_proteins_per_genome"] = sum(protein_counts) / len(protein_counts)
    
    stats["sample_size"] = sample_size
    return stats

def get_prokka_version():
    """Get Prokka version information."""
    try:
        import subprocess
        # First try to activate conda environment and get version
        conda_cmd = "eval \"$(/home/pmuench/miniconda3/bin/conda shell.bash hook)\" && conda activate efs_diversity && prokka --version 2>&1"
        result = subprocess.run(conda_cmd, shell=True, 
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0 and result.stdout.strip():
            # Prokka outputs version to stderr, but we're capturing both with 2>&1
            version_line = result.stdout.strip().split('\n')[-1]
            if 'prokka' in version_line.lower():
                return version_line
        
        # Fallback to direct prokka call
        result = subprocess.run(['prokka', '--version'], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            # Check both stdout and stderr
            output = result.stdout.strip() if result.stdout.strip() else result.stderr.strip()
            if output:
                return output.split('\n')[0]
    except:
        pass
    return "Prokka version not available"

def generate_manuscript_stats(output_file=None):
    """Generate all statistics for manuscript.
    
    Args:
        output_file (str): Optional path to save output to file. If None, prints to console.
    """
    
    # Collect all output lines
    output_lines = []
    
    def add_line(text=""):
        """Add a line to both console and file output."""
        print(text)
        output_lines.append(text)
    
    add_line("Prokka Annotation Statistics for Manuscript")
    add_line("=" * 50)
    
    # Input genomes
    input_count = count_input_genomes()
    add_line(f"\n1. Input Genomes:")
    add_line(f"   Total E. faecalis assemblies: {input_count:,}")
    
    # Prokka version
    version = get_prokka_version()
    add_line(f"\n2. Software Version:")
    add_line(f"   {version}")
    
    # Prokka outputs
    add_line(f"\n3. Annotation Results:")
    prokka_stats = analyze_prokka_outputs()
    
    if "error" not in prokka_stats:
        add_line(f"   Total genomes processed: {prokka_stats['total_genomes']:,}")
        add_line(f"   Sample analyzed: {prokka_stats.get('sample_size', 0):,}")
        add_line(f"   Complete annotations: {prokka_stats['complete_genomes']:,}")
        if prokka_stats['complete_genomes'] > 0:
            completion_rate = prokka_stats['complete_genomes'] / prokka_stats.get('sample_size', 1) * 100
            add_line(f"   Success rate: {completion_rate:.1f}%")
        
        if prokka_stats['avg_genes_per_genome'] > 0:
            add_line(f"   Average genes per genome: {prokka_stats['avg_genes_per_genome']:.0f}")
            add_line(f"   Median genes per genome: {prokka_stats.get('median_genes_per_genome', 0):.0f}")
            add_line(f"   Gene count range: {prokka_stats.get('min_genes_per_genome', 0):.0f}-{prokka_stats.get('max_genes_per_genome', 0):.0f}")
            add_line(f"   Standard deviation: {prokka_stats.get('std_genes_per_genome', 0):.0f}")
            add_line(f"   Total genes annotated: {prokka_stats['total_genes']:,}")
            
            # Outlier analysis
            if prokka_stats.get('num_very_low_gene', 0) > 0:
                add_line(f"\n   Quality Control - Potential Outliers:")
                add_line(f"   Genomes with <1000 genes: {prokka_stats['num_very_low_gene']:,}")
                add_line(f"   Statistical outliers: {prokka_stats.get('num_outliers', 0):,}")
                
                if prokka_stats.get('very_low_gene_genomes'):
                    add_line(f"   Examples of low-gene genomes:")
                    for genome_id, gene_count in prokka_stats['very_low_gene_genomes'][:5]:
                        add_line(f"     {genome_id}: {gene_count} genes")
                
                # Percentile distribution
                if prokka_stats.get('gene_count_percentiles'):
                    add_line(f"   Gene count distribution:")
                    percentiles = prokka_stats['gene_count_percentiles']
                    add_line(f"     1%: {percentiles['1%']:.0f}, 5%: {percentiles['5%']:.0f}, 10%: {percentiles['10%']:.0f}")
                    add_line(f"     25%: {percentiles['25%']:.0f}, 50%: {percentiles['50%']:.0f}, 75%: {percentiles['75%']:.0f}")
                    add_line(f"     90%: {percentiles['90%']:.0f}, 95%: {percentiles['95%']:.0f}, 99%: {percentiles['99%']:.0f}")
    else:
        add_line(f"   {prokka_stats['error']}")
    
    # Output file types
    add_line(f"\n4. Output Files Generated per Genome:")
    file_types = {
        "gff": "Gene feature format (annotations)",
        "faa": "Protein sequences (amino acids)", 
        "ffn": "Nucleotide sequences (genes)",
        "fna": "Nucleotide sequences (contigs)",
        "gbk": "GenBank format",
        "tsv": "Feature table",
        "txt": "Statistics summary",
        "log": "Processing log"
    }
    
    for ext, description in file_types.items():
        add_line(f"   .{ext}: {description}")
    
    # Summary for manuscript
    add_line("\n" + "=" * 50)
    add_line("MANUSCRIPT NUMBERS SUMMARY:")
    add_line("=" * 50)
    add_line(f"Input assemblies: {input_count:,} E. faecalis genomes")
    if "error" not in prokka_stats and prokka_stats['avg_genes_per_genome'] > 0:
        add_line(f"Successfully annotated: {prokka_stats['complete_genomes']:,} genomes")
        add_line(f"Average genes per genome: {prokka_stats['avg_genes_per_genome']:.0f}")
        add_line(f"Total genes predicted: {prokka_stats['total_genes']:,}")
        add_line(f"Output file types: 8 per genome (.gff, .faa, .ffn, .fna, .gbk, .tsv, .txt, .log)")
    
    # Write to file if specified
    if output_file:
        try:
            with open(output_file, 'w') as f:
                for line in output_lines:
                    f.write(line + '\n')
            add_line(f"\nResults saved to: {output_file}")
        except Exception as e:
            add_line(f"\nError saving to file: {e}")
    
    return {
        'input_count': input_count,
        'prokka_stats': prokka_stats,
        'version': version,
        'output_lines': output_lines
    }

if __name__ == "__main__":
    import sys
    
    # Change to script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    
    # Check for output file argument
    output_file = None
    if len(sys.argv) > 1:
        output_file = sys.argv[1]
    
    stats = generate_manuscript_stats(output_file)
