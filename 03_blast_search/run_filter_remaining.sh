#!/bin/bash
#SBATCH --job-name=filter_remaining
#SBATCH --output=logs/filter_remaining_%j.out
#SBATCH --error=logs/filter_remaining_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu

echo "=============================================="
echo "Filtering Remaining Output Files"
echo "Started: $(date)"
echo "=============================================="

# Create logs directory
mkdir -p logs

# Activate conda
eval "$(/home/pmuench/miniconda3/bin/conda shell.bash hook)"
conda activate efs_diversity

cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/03_blast_search

# Run the filter script
python filter_remaining_outputs.py

echo ""
echo "=============================================="
echo "Completed: $(date)"
echo "=============================================="
