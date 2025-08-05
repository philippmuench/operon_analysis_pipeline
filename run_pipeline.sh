#!/bin/bash
#SBATCH --job-name=operon_pipeline
#SBATCH --output=pipeline_%j.out
#SBATCH --error=pipeline_%j.err
#SBATCH --time=72:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu

# Master script to run the complete operon analysis pipeline
# Each step submits its own jobs as needed

set -e  # Exit on error

echo "========================================="
echo "E. faecalis Operon Analysis Pipeline"
echo "Started: $(date)"
echo "========================================="

# Load environment
source /vol/projects/BIFO/utils/loadEnv
conda activate efs_diversity

# Base directory
BASE_DIR=$(pwd)

# Step 1: Prokka annotation
echo -e "\n[Step 1] Running Prokka annotation..."
cd 01_prokka_annotation

if [ $(find ../prokka_output -name "*.gff" 2>/dev/null | wc -l) -lt 8000 ]; then
    echo "Submitting Prokka jobs..."
    PROKKA_JOB=$(sbatch --parsable run_prokka.sh)
    echo "Prokka job: $PROKKA_JOB"
    
    # Wait for completion
    while squeue -j $PROKKA_JOB 2>/dev/null | grep -q $PROKKA_JOB; do
        sleep 300  # Check every 5 minutes
        echo "  Prokka progress: $(./check_progress.sh | grep Progress)"
    done
else
    echo "Prokka annotation already complete"
fi

cd $BASE_DIR

# Step 2: Extract operon genes
echo -e "\n[Step 2] Extracting operon genes..."
cd 02_operon_extraction

if [ ! -f ../data/operon_genes.fasta ]; then
    python extract_operon_genes.py
else
    echo "Operon genes already extracted"
fi

cd $BASE_DIR

# Step 3: BLAST search
echo -e "\n[Step 3] Running BLAST search..."
cd 03_blast_search

if [ ! -f ../results/operon_presence_summary.csv ]; then
    echo "Submitting BLAST jobs..."
    BLAST_JOB=$(sbatch --parsable run_blast_search.sh)
    echo "BLAST job: $BLAST_JOB"
    
    # Wait for completion
    while squeue -j $BLAST_JOB 2>/dev/null | grep -q $BLAST_JOB; do
        sleep 300
    done
    
    echo "Processing BLAST results..."
    python process_blast_results.py
else
    echo "BLAST search already complete"
fi

cd $BASE_DIR

# Step 4: Core gene analysis
echo -e "\n[Step 4] Analyzing core genes..."
cd 04_core_gene_analysis

if [ ! -f ../data/core_genes_95pct.txt ]; then
    python identify_core_genes.py
fi

if [ ! -f ../results/core_gene_diversity.csv ]; then
    echo "Submitting core gene analysis..."
    CORE_JOB=$(sbatch --parsable run_core_analysis.sh)
    echo "Core analysis job: $CORE_JOB"
    
    while squeue -j $CORE_JOB 2>/dev/null | grep -q $CORE_JOB; do
        sleep 300
    done
fi

cd $BASE_DIR

# Step 5: Diversity analysis
echo -e "\n[Step 5] Running diversity analysis..."
cd 05_diversity_analysis

if [ ! -f ../results/substitution_analysis.csv ]; then
    echo "Submitting diversity analysis..."
    DIV_JOB=$(sbatch --parsable run_diversity_analysis.sh)
    echo "Diversity job: $DIV_JOB"
    
    while squeue -j $DIV_JOB 2>/dev/null | grep -q $DIV_JOB; do
        sleep 300
    done
fi

cd $BASE_DIR

# Step 6: Generate visualizations
echo -e "\n[Step 6] Generating visualizations..."
cd 06_visualization

python create_all_figures.py

cd $BASE_DIR

echo -e "\n========================================="
echo "Pipeline completed: $(date)"
echo "========================================="
echo "Results available in:"
echo "  - results/operon_presence_summary.csv"
echo "  - results/core_gene_diversity.csv"
echo "  - results/figures/"
echo "========================================="