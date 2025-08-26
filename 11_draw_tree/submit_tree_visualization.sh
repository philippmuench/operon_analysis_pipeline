#!/bin/bash
#SBATCH --job-name=draw_tree
#SBATCH --output=draw_tree_%j.log
#SBATCH --error=draw_tree_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --partition=cpu

echo "Starting tree visualization job"
echo "Job ID: $SLURM_JOB_ID"
echo "Date: $(date)"
echo "================================"

# Activate conda environment
eval "$(/home/pmuench/miniconda3/bin/conda shell.bash hook)"
conda activate efs_diversity

# Run the visualization pipeline
cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/11_draw_tree
python3 draw_tree_with_metadata.py

echo "================================"
echo "Job completed at: $(date)"