#!/bin/bash
#SBATCH --job-name=operon_variation_slices
#SBATCH --output=operon_variation_slices_%j.out
#SBATCH --error=operon_variation_slices_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=cpu

set -euo pipefail

cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/13e_full_operon_mapping

eval "$(/home/pmuench/miniconda3/bin/conda shell.bash hook)"
conda activate efs_diversity

python plot_variation_slices.py
