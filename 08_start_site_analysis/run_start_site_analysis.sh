#!/bin/bash
#SBATCH --job-name=start_site_analysis
#SBATCH --output=start_site_%j.out
#SBATCH --error=start_site_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8
#SBATCH --partition=cpu

set -euo pipefail

# Try to load common env if present, but don't fail if missing
if [ -f /vol/projects/BIFO/utils/loadEnv ]; then
  source /vol/projects/BIFO/utils/loadEnv || true
else
  echo "Warning: /vol/projects/BIFO/utils/loadEnv not found; continuing without it" >&2
fi

# Try to activate conda env if available
if command -v conda >/dev/null 2>&1; then
  # Initialize conda for non-interactive shells if needed
  if [ -z "${CONDA_EXE:-}" ]; then
    CONDA_BASE=$(conda info --base 2>/dev/null || true)
    if [ -n "$CONDA_BASE" ] && [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
      # shellcheck disable=SC1090
      source "$CONDA_BASE/etc/profile.d/conda.sh" || true
    fi
  fi
  conda activate efs_diversity 2>/dev/null || echo "Warning: could not activate conda env 'efs_diversity'; using system python" >&2
else
  echo "Warning: conda not found in PATH; using system python" >&2
fi

PROKKA_DIR=${1:-"../01_prokka_annotation/output/prokka_results"}
GENE_FASTA=${2:-"../02_reference_operon_extraction/output/operon_genes_nt.fasta"}
OUT_DIR=${3:-"./output"}

mkdir -p "$OUT_DIR"

echo "Running start-site analysis on $PROKKA_DIR"
python analyze_start_sites.py \
  --prokka_dir "$PROKKA_DIR" \
  --gene_reference_fasta "$GENE_FASTA" \
  --output_dir "$OUT_DIR" \
  --blast_dir "../03_blast_search/output/blast_results_prokka_variants" \
  --max_workers "$SLURM_CPUS_PER_TASK"

echo "Done. Output in $OUT_DIR"


