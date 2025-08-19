#!/bin/bash
#SBATCH --job-name=ptsA_upstream
#SBATCH --output=ptsA_upstream_%j.out
#SBATCH --error=ptsA_upstream_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8
#SBATCH --partition=cpu

set -euo pipefail

# Optional environment load
if [ -f /vol/projects/BIFO/utils/loadEnv ]; then
  source /vol/projects/BIFO/utils/loadEnv || true
else
  echo "Warning: /vol/projects/BIFO/utils/loadEnv not found; continuing without it" >&2
fi

# Try to activate conda env if available
if command -v conda >/dev/null 2>&1; then
  if [ -z "${CONDA_EXE:-}" ]; then
    CONDA_BASE=$(conda info --base 2>/dev/null || true)
    if [ -n "$CONDA_BASE" ] && [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
      source "$CONDA_BASE/etc/profile.d/conda.sh" || true
    fi
  fi
  conda activate efs_diversity 2>/dev/null || echo "Warning: could not activate conda env 'efs_diversity'; using system python" >&2
else
  echo "Warning: conda not found in PATH; using system python" >&2
fi

PROKKA_DIR=${1:-"../01_prokka_annotation/output/prokka_results"}
BLAST_DIR=${2:-"../03_blast_search/output/blast_results_prokka_variants"}
SUMMARY_TSV=${3:-"./output/start_site_summary.tsv"}
OUT_DIR=${4:-"./output/ptsA_upstream"}
# Optional toggles/params
FULL_GENE=${5:-"false"}          # "true" to include entire CDS downstream of start in windows
UPSTREAM_LEN=${6:-120}
DOWNSTREAM_LEN=${7:-60}

mkdir -p "$OUT_DIR"

# Ensure MAFFT_TMPDIR is set and usable (mirror policy from other modules)
# Prefer an explicit cluster scratch path; allow user override
export MAFFT_TMPDIR="${MAFFT_TMPDIR:-/vol/tmp}"
if [ -z "$MAFFT_TMPDIR" ]; then
  echo "ERROR: MAFFT_TMPDIR is not set. Please export MAFFT_TMPDIR (e.g., /vol/tmp) before running." >&2
  exit 1
fi
if [ ! -d "$MAFFT_TMPDIR" ] || [ ! -w "$MAFFT_TMPDIR" ]; then
  echo "ERROR: MAFFT_TMPDIR=$MAFFT_TMPDIR is not an existing writable directory." >&2
  exit 1
fi
echo "Using MAFFT_TMPDIR=$MAFFT_TMPDIR"

echo "Running ptsA upstream analysis"
echo "  full_gene_downstream=${FULL_GENE} | upstream_len=${UPSTREAM_LEN} | downstream_len=${DOWNSTREAM_LEN}"
# Build optional flag strings
ANALYZE_FLAGS=()
PLOT_FLAGS=()
if [ "${FULL_GENE}" = "true" ]; then
  ANALYZE_FLAGS+=("--full_gene_downstream")
  PLOT_FLAGS+=("--full-gene")
fi

python analyze_ptsA_upstream.py \
  --prokka_dir "$PROKKA_DIR" \
  --blast_dir "$BLAST_DIR" \
  --summary_tsv "$SUMMARY_TSV" \
  --out_dir "$OUT_DIR" \
  --upstream_len "$UPSTREAM_LEN" \
  --downstream_len "$DOWNSTREAM_LEN" \
  --threads "${SLURM_CPUS_PER_TASK:-8}" \
  "${ANALYZE_FLAGS[@]}"

echo "Done. See $OUT_DIR"

# If alignment exists, generate heatmap and profiles with start-site overlay
MSA="$OUT_DIR/ptsA_windows.aln.fa"
RAW="$OUT_DIR/ptsA_windows.fa"
if [ -s "$MSA" ] && [ -s "$RAW" ]; then
  echo "Generating alignment heatmap and profiles"
  python plot_alignment_profiles.py \
    --msa "$MSA" \
    --outdir "$OUT_DIR" \
    --raw-windows "$RAW" \
    --upstream-len "$UPSTREAM_LEN" \
    --downstream-len "$DOWNSTREAM_LEN" \
    --max-rows 200 \
    --save-pdf \
    "${PLOT_FLAGS[@]}"
fi


