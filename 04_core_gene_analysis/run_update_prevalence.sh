#!/bin/bash
#SBATCH --job-name=update_prevalence
#SBATCH --output=logs/update_prevalence_%j.out
#SBATCH --error=logs/update_prevalence_%j.err
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu

echo "=============================================="
echo "Updating Prevalence Stats and Threshold Curve"
echo "Started: $(date)"
echo "=============================================="

mkdir -p logs

eval "$(/home/pmuench/miniconda3/bin/conda shell.bash hook)"
conda activate efs_diversity

cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/04_core_gene_analysis

python update_prevalence_from_msas.py

echo ""
echo "=============================================="
echo "Completed: $(date)"
echo "=============================================="
