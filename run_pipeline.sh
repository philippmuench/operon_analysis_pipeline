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
cd 02_reference_operon_extraction

python extract_operon_genes.py

cd $BASE_DIR

# Step 3: BLAST search
echo -e "\n[Step 3] Running BLAST search..."
cd 03_blast_search

echo "Submitting BLAST jobs..."
BLAST_JOB=$(sbatch --parsable run_blast_search.sh)
echo "BLAST job: $BLAST_JOB"

# Wait for completion
while squeue -j $BLAST_JOB 2>/dev/null | grep -q $BLAST_JOB; do
    sleep 300
done

echo "Processing BLAST results..."
python process_blast_results.py

cd $BASE_DIR

# Step 4: Core gene analysis
echo -e "\n[Step 4] Analyzing core genes..."
cd 04_core_gene_analysis

echo "Submitting core gene analysis..."
CORE_JOB=$(sbatch --parsable run_core_analysis.sh)
echo "Core analysis job: $CORE_JOB"

while squeue -j $CORE_JOB 2>/dev/null | grep -q $CORE_JOB; do
    sleep 300
done

cd $BASE_DIR

echo -e "\n[Step 5] Operon assembly-based extraction and MSAs..."
cd 05_operon_assembly_extraction
OP_JOB=$(sbatch --parsable run_operon_extraction.sh)
echo "Operon extraction job: $OP_JOB"
while squeue -j $OP_JOB 2>/dev/null | grep -q $OP_JOB; do
    sleep 300
done

cd $BASE_DIR

echo -e "\n[Step 6] Diversity analysis..."
cd 06_diversity_analysis
DIV_JOB=$(sbatch --parsable run_complete_diversity_analysis.sh)
echo "Diversity job: $DIV_JOB"
while squeue -j $DIV_JOB 2>/dev/null | grep -q $DIV_JOB; do
    sleep 300
done

cd $BASE_DIR

echo -e "\n[Step 7] dN/dS analysis..."
cd 07_dnds_analysis
DNDS_JOB=$(sbatch --parsable run_dnds_analysis.sh)
echo "dN/dS job: $DNDS_JOB"
while squeue -j $DNDS_JOB 2>/dev/null | grep -q $DNDS_JOB; do
    sleep 300
done

echo -e "\n========================================="
echo "Pipeline completed: $(date)"
echo "========================================="
echo "Results available in:"
echo "  - 03_blast_search/output/operon_presence_summary.csv"
echo "  - 04_core_gene_analysis/output/gene_prevalence_stats.csv"
echo "  - 05_operon_assembly_extraction/output/"
echo "  - 06_diversity_analysis/output/"
echo "  - 07_dnds_analysis/output/"
echo "========================================="