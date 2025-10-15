#!/bin/bash
#SBATCH --job-name=codon_variation_viz
#SBATCH --output=codon_variation_viz_%j.out
#SBATCH --error=codon_variation_viz_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --partition=cpu

set -euo pipefail

ALIGN_DIR="../05_operon_assembly_extraction/output/msa/dna_alignments"
OUTPUT_DIR="output"
PY_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --alignment-dir)
            ALIGN_DIR="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        *)
            PY_ARGS+=("$1")
            shift
            ;;
    esac
done

echo "=========================================="
echo "Codon Variation Visualisation"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "Working directory: $(pwd)"
echo "Script args: $*"
echo "=========================================="

# Activate conda environment
echo ""
echo "Setting up environment..."
eval "$(${HOME}/miniconda3/bin/conda shell.bash hook)"
conda activate efs_diversity

echo "Active environment: $CONDA_DEFAULT_ENV"
python --version

if [ ! -d "$ALIGN_DIR" ]; then
    echo "✗ Alignment directory not found: $ALIGN_DIR"
    exit 1
fi

echo ""
echo "Output directory: $(pwd)/$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"/{tables,plots}

echo ""
echo "Running visualisation..."
python visualize_codon_variation.py \
    --alignment-dir "$ALIGN_DIR" \
    --output-dir "$OUTPUT_DIR" \
    "${PY_ARGS[@]}"

EXIT_STATUS=$?
if [ $EXIT_STATUS -ne 0 ]; then
    echo ""
    echo "=========================================="
    echo "ERROR: Visualisation failed (exit $EXIT_STATUS)"
    echo "=========================================="
    exit $EXIT_STATUS
fi

echo ""
echo "=========================================="
echo "Visualisation completed successfully"
echo "=========================================="

echo "Generated tables:"
ls -1 "$OUTPUT_DIR"/tables 2>/dev/null || echo "  (none)"

echo "Generated plots:"
ls -1 "$OUTPUT_DIR"/plots 2>/dev/null || echo "  (none)"

echo "Job finished at $(date)"
