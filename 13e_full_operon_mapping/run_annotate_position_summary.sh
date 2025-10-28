#!/bin/bash
#SBATCH --job-name=operon_annotate_positions
#SBATCH --output=operon_annotate_positions_%j.out
#SBATCH --error=operon_annotate_positions_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=cpu

set -euo pipefail

cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/13e_full_operon_mapping

eval "$(/home/pmuench/miniconda3/bin/conda shell.bash hook)"
conda activate efs_diversity

python annotate_position_summary.py \
    --combined output_full_run/position_stats/combined_position_summary.tsv \
    --annotation-table ../13a_new_reference_operon_extraction/output/operon_genes.tsv \
    --variation-table ../13a_new_reference_operon_extraction/output/operon_variations.tsv
