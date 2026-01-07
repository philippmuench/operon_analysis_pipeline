#!/bin/bash
#SBATCH --job-name=regen_metrics
#SBATCH --output=logs/regen_metrics_%j.out
#SBATCH --error=logs/regen_metrics_%j.err
#SBATCH --time=00:15:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu

echo "=============================================="
echo "Regenerating Conservation Metrics and Plots"
echo "Started: $(date)"
echo "=============================================="

mkdir -p logs

eval "$(/home/pmuench/miniconda3/bin/conda shell.bash hook)"
conda activate efs_diversity

cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/05_operon_assembly_extraction

python regenerate_metrics_plots.py

echo ""
echo "=============================================="
echo "Completed: $(date)"
echo "=============================================="
