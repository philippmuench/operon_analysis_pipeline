#!/bin/bash
#SBATCH --job-name=diversity_all
#SBATCH --output=diversity_%j.out
#SBATCH --error=diversity_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=32
#SBATCH --partition=cpu

# Complete diversity analysis pipeline

source /vol/projects/BIFO/utils/loadEnv
conda activate efs_diversity

echo "=== Diversity Analysis Pipeline ==="
echo "Started: $(date)"

# Paths
MSA_DIR="../../msa_output"
CORE_SEQS="../../core_gene_sequences_95pct"
OUTPUT_DIR="../../results"
mkdir -p $OUTPUT_DIR

# Step 1: Calculate core gene diversity
echo -e "\n[Step 1] Calculating core gene diversity..."
if [ ! -f "$OUTPUT_DIR/core_gene_diversity_all.csv" ]; then
    python calculate_core_gene_diversity.py \
        --input_dir $CORE_SEQS \
        --output $OUTPUT_DIR/core_gene_diversity_all.csv \
        --cores 32
else
    echo "Core gene diversity already calculated"
fi

# Step 2: Analyze operon gene diversity
echo -e "\n[Step 2] Analyzing operon gene diversity..."
if [ ! -f "$OUTPUT_DIR/operon_diversity_detailed.csv" ]; then
    python analyze_msa_diversity.py \
        --msa_dir $MSA_DIR \
        --output $OUTPUT_DIR/operon_diversity_detailed.csv
else
    echo "Operon diversity already calculated"
fi

# Step 3: Calculate dN/dS for operon genes
echo -e "\n[Step 3] Calculating dN/dS for operon genes..."
if [ ! -f "$OUTPUT_DIR/operon_dnds_results.csv" ]; then
    python analyze_operon_substitutions.py \
        --msa_dir $MSA_DIR \
        --output $OUTPUT_DIR/operon_dnds_results.csv
else
    echo "Operon dN/dS already calculated"
fi

# Step 4: Compare operon vs core genes (dN/dS)
echo -e "\n[Step 4] Comparing substitution patterns..."
if [ ! -f "$OUTPUT_DIR/operon_vs_core_dnds.csv" ]; then
    python compare_substitutions_all.py \
        --operon_dir $MSA_DIR \
        --core_dir $CORE_SEQS \
        --output $OUTPUT_DIR/operon_vs_core_dnds.csv \
        --cores 32
else
    echo "Substitution comparison already complete"
fi

# Step 5: Generate summary statistics
echo -e "\n[Step 5] Generating summary statistics..."
python generate_diversity_summary.py \
    --core_diversity $OUTPUT_DIR/core_gene_diversity_all.csv \
    --operon_diversity $OUTPUT_DIR/operon_diversity_detailed.csv \
    --dnds_results $OUTPUT_DIR/operon_vs_core_dnds.csv \
    --output $OUTPUT_DIR/diversity_summary.json

echo -e "\nDiversity analysis completed: $(date)"
echo "Results in: $OUTPUT_DIR"