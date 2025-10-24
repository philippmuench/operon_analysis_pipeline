#!/bin/bash
#SBATCH --job-name=operon_raw_stats
#SBATCH --output=operon_raw_stats_%j.out
#SBATCH --error=operon_raw_stats_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=cpu

set -euo pipefail

cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/13e_full_operon_mapping

eval "$(/home/pmuench/miniconda3/bin/conda shell.bash hook)"
conda activate efs_diversity

python process_raw_blast.py \
    --legacy-raw output_full_run/raw_blast/legacy \
    --updated-raw output_full_run/raw_blast/updated \
    --legacy-reference ../02_reference_operon_extraction/operon.gb \
    --updated-reference ../13a_new_reference_operon_extraction/FL_operon_with_SNPs.gb \
    --min-identity 90 \
    --min-coverage 80 \
    --output output_full_run/position_stats
