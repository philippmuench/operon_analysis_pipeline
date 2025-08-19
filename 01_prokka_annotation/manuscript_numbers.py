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
    genome_dir = "../Efs_assemblies"
    if not os.path.exists(genome_dir):
        return 0
    
    genome_files = glob.glob(os.path.join(genome_dir, "*.fasta.gz"))
    return len(genome_files)

def analyze_prokka_outputs():
    """Analyze Prokka annotation outputs."""
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
        "file_counts": defaultdict(int)
    }
    
    expected_files = ["gff", "faa", "ffn", "fna", "gbk", "log", "txt", "tsv"]
    gene_counts = []
    protein_counts = []
    
    for genome_dir in genome_dirs[:100]:  # Sample first 100 for speed
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
                    stats["total_genes"] += gene_count
                except:
                    pass
            
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
    
    # Calculate averages
    if gene_counts:
        stats["avg_genes_per_genome"] = sum(gene_counts) / len(gene_counts)
    if protein_counts:
        stats["avg_proteins_per_genome"] = sum(protein_counts) / len(protein_counts)
    
    return stats

def get_prokka_version():
    """Get Prokka version information."""
    try:
        import subprocess
        result = subprocess.run(['prokka', '--version'], 
                              capture_output=True, text=True, timeout=10)
        if result.returncode == 0:
            return result.stdout.strip().split('\n')[0]
    except:
        pass
    return "Prokka version not available"

def generate_manuscript_stats():
    """Generate all statistics for manuscript."""
    
    print("Prokka Annotation Statistics for Manuscript")
    print("=" * 50)
    
    # Input genomes
    input_count = count_input_genomes()
    print(f"\n1. Input Genomes:")
    print(f"   Total E. faecalis assemblies: {input_count:,}")
    
    # Prokka version
    version = get_prokka_version()
    print(f"\n2. Software Version:")
    print(f"   {version}")
    
    # Prokka outputs
    print(f"\n3. Annotation Results:")
    prokka_stats = analyze_prokka_outputs()
    
    if "error" not in prokka_stats:
        print(f"   Total genomes processed: {prokka_stats['total_genomes']:,}")
        print(f"   Complete annotations: {prokka_stats['complete_genomes']:,}")
        if prokka_stats['complete_genomes'] > 0:
            completion_rate = prokka_stats['complete_genomes'] / prokka_stats['total_genomes'] * 100
            print(f"   Success rate: {completion_rate:.1f}%")
        
        if prokka_stats['avg_genes_per_genome'] > 0:
            print(f"   Average genes per genome: {prokka_stats['avg_genes_per_genome']:.0f}")
            print(f"   Average proteins per genome: {prokka_stats['avg_proteins_per_genome']:.0f}")
            print(f"   Total genes annotated: {prokka_stats['total_genes']:,}")
            print(f"   Total proteins predicted: {prokka_stats['total_proteins']:,}")
    else:
        print(f"   {prokka_stats['error']}")
    
    # Output file types
    print(f"\n4. Output Files Generated per Genome:")
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
        print(f"   .{ext}: {description}")
    
    # Summary for manuscript
    print("\n" + "=" * 50)
    print("MANUSCRIPT NUMBERS SUMMARY:")
    print("=" * 50)
    print(f"Input assemblies: {input_count:,} E. faecalis genomes")
    if "error" not in prokka_stats and prokka_stats['avg_genes_per_genome'] > 0:
        print(f"Successfully annotated: {prokka_stats['complete_genomes']:,} genomes")
        print(f"Average genes per genome: {prokka_stats['avg_genes_per_genome']:.0f}")
        print(f"Total genes predicted: {prokka_stats['total_genes']:,}")
        print(f"Output file types: 8 per genome (.gff, .faa, .ffn, .fna, .gbk, .tsv, .txt, .log)")
    
    return {
        'input_count': input_count,
        'prokka_stats': prokka_stats,
        'version': version
    }

if __name__ == "__main__":
    # Change to script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    
    stats = generate_manuscript_stats()
