#!/bin/bash
#SBATCH --job-name=dnds_analysis
#SBATCH --output=dnds_%j.out
#SBATCH --error=dnds_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpu

# Run dN/dS analysis for operon and core genes

source /vol/projects/BIFO/utils/loadEnv
conda activate efs_diversity

echo "=== dN/dS Analysis ==="
echo "Started: $(date)"

# Create output directory
mkdir -p ../../results/dnds_analysis

# Step 1: Calculate dN/dS for operon genes
echo -e "\n[1] Analyzing operon genes..."
python analyze_dnds.py \
    --input_dir ../../msa_output \
    --output ../../results/dnds_analysis/operon_dnds.csv \
    --cores 16

# Step 2: Calculate dN/dS for core genes (sample if too many)
echo -e "\n[2] Analyzing core genes..."
# First check how many alignments we have
n_alignments=$(ls ../../core_gene_sequences_95pct/*.fasta 2>/dev/null | wc -l)

if [ $n_alignments -gt 500 ]; then
    echo "Sampling 500 core genes for dN/dS analysis..."
    # Create temp directory with sampled genes
    mkdir -p temp_core_sample
    ls ../../core_gene_sequences_95pct/*.fasta | shuf -n 500 | while read f; do
        cp "$f" temp_core_sample/
    done
    
    python analyze_dnds.py \
        --input_dir temp_core_sample \
        --output ../../results/dnds_analysis/core_genes_dnds.csv \
        --cores 16
    
    rm -rf temp_core_sample
else
    python analyze_dnds.py \
        --input_dir ../../core_gene_sequences_95pct \
        --output ../../results/dnds_analysis/core_genes_dnds.csv \
        --cores 16
fi

# Step 3: Compare results
echo -e "\n[3] Comparing operon vs core genes..."
python compare_dnds_results.py \
    --operon_file ../../results/dnds_analysis/operon_dnds.csv \
    --core_file ../../results/dnds_analysis/core_genes_dnds.csv \
    --output ../../results/dnds_analysis/comparison_summary.txt

echo -e "\nAnalysis completed: $(date)"
echo "Results in: ../../results/dnds_analysis/"