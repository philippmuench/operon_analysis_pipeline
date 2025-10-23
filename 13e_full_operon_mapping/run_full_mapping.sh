#!/bin/bash
#SBATCH --job-name=full_operon_mapping
#SBATCH --output=full_operon_mapping_%j.out
#SBATCH --error=full_operon_mapping_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=cpu

set -euo pipefail

cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/13e_full_operon_mapping

eval "$(/home/pmuench/miniconda3/bin/conda shell.bash hook)"
conda activate efs_diversity

python compare_full_operon_mapping.py \
    --legacy-gb ../02_reference_operon_extraction/operon.gb \
    --updated-gb ../13a_new_reference_operon_extraction/FL_operon_with_SNPs.gb \
    --assemblies-dir ../../Efs_assemblies \
    --threads "$SLURM_CPUS_PER_TASK" \
    --output output_full_run \
    --min-coverage 80 \
    --min-identity 90

