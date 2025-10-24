#!/bin/bash
#SBATCH --job-name=operon_variant_summary
#SBATCH --output=operon_variant_summary_%j.out
#SBATCH --error=operon_variant_summary_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=cpu

set -euo pipefail

cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis

eval "$(/home/pmuench/miniconda3/bin/conda shell.bash hook)"
conda activate efs_diversity

./13e_full_operon_mapping/analyze_updated_operon_variants.py \
  --site-summary 13e_full_operon_mapping/output_full_run/variant_site_summary.tsv \
  --site-table 13e_full_operon_mapping/output_full_run/variant_site_allele_counts.tsv \
  --site-calls 13e_full_operon_mapping/output_full_run/variant_site_calls.tsv \
  --figure-out 13e_full_operon_mapping/output_full_run/variant_site_counts.png
