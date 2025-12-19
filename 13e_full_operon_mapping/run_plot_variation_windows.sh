#!/bin/bash
#SBATCH --job-name=operon_variation_windows
#SBATCH --output=operon_variation_windows_%j.out
#SBATCH --error=operon_variation_windows_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=cpu

set -euo pipefail

cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/13e_full_operon_mapping

eval "$(/home/pmuench/miniconda3/bin/conda shell.bash hook)"
conda activate efs_diversity

python plot_variation_windows.py \
    --combined output_full_run/position_stats/combined_position_summary.tsv \
    --window 10 \
    --output-dir output_full_run/position_stats/variation_windows \
    --legacy-label Legacy \
    --updated-label Updated
