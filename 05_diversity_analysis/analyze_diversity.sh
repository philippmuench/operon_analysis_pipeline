#!/bin/bash
#SBATCH --job-name=diversity
#SBATCH --output=logs/diversity_%j.out
#SBATCH --error=logs/diversity_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpu

# Wait for extraction array job
if [ -f .extraction_job_id ]; then
    EXTRACT_JOB=$(cat .extraction_job_id)
    echo "Waiting for extraction job $EXTRACT_JOB to complete..."
    while squeue -j $EXTRACT_JOB 2>/dev/null | grep -q $EXTRACT_JOB; do
        sleep 30
    done
fi

echo "Starting diversity analysis at $(date)"

# Check extraction results
N_GENES=$(ls -1 core_gene_sequences/*.fasta 2>/dev/null | wc -l)
echo "Found $N_GENES gene sequence files"

# Calculate diversity
python calculate_core_gene_diversity.py

# Create visualizations
if [ $? -eq 0 ]; then
    python visualize_diversity_comparison.py
    echo "Analysis completed successfully at $(date)"
else
    echo "Error in diversity calculation"
    exit 1
fi
