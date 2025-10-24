#!/bin/bash
#SBATCH --job-name=operon_position_plots
#SBATCH --output=operon_position_plots_%j.out
#SBATCH --error=operon_position_plots_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=cpu

set -euo pipefail

cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/13e_full_operon_mapping

eval "$(/home/pmuench/miniconda3/bin/conda shell.bash hook)"
conda activate efs_diversity

python generate_position_plots.py \
    --input output_full_run/position_stats/combined_position_summary.tsv \
    --output-dir output_full_run/position_stats/plots \
    --legacy-annotation ../02_reference_operon_extraction/output/operon_genes.tsv \
    --updated-annotation ../13a_new_reference_operon_extraction/output/operon_genes.tsv
