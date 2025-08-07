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
mkdir -p output

# Step 1: Calculate dN/dS for operon genes
echo -e "\n[1] Analyzing operon genes..."
python analyze_dnds.py \
    --input_dir ../05_operon_assembly_extraction/output/msa/dna_alignments \
    --output output/operon_dnds.csv \
    --cores 16 \
    --file_pattern "*_aligned.fasta"

# Step 2: Calculate dN/dS for core genes (sample if too many)
echo -e "\n[2] Analyzing core genes..."
# First check how many alignments we have
n_alignments=$(ls ../04_core_gene_analysis/output/core_gene_alignments/*_aligned.fasta 2>/dev/null | wc -l)

if [ $n_alignments -gt 500 ]; then
    echo "Sampling 500 core genes for dN/dS analysis..."
    # Create temp directory with sampled genes
    mkdir -p temp_core_sample
    ls ../04_core_gene_analysis/output/core_gene_alignments/*_aligned.fasta | shuf -n 500 | while read f; do
        cp "$f" temp_core_sample/
    done
    
    python analyze_dnds.py \
        --input_dir temp_core_sample \
        --output output/core_genes_dnds.csv \
        --cores 16 \
        --file_pattern "*_aligned.fasta"
    
    rm -rf temp_core_sample
else
    python analyze_dnds.py \
        --input_dir ../04_core_gene_analysis/output/core_gene_alignments \
        --output output/core_genes_dnds.csv \
        --cores 16 \
        --file_pattern "*_aligned.fasta"
fi

# Step 3: Compare results
echo -e "\n[3] Comparing operon vs core genes..."
python compare_dnds_results.py \
    --operon_file output/operon_dnds.csv \
    --core_file output/core_genes_dnds.csv \
    --output output/comparison_summary.txt

echo -e "\nAnalysis completed: $(date)"
echo "Results in: $(pwd)/output/"