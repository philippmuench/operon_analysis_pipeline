#!/bin/bash
#SBATCH --job-name=blast_pipeline
#SBATCH --output=output/pipeline_%A_%a.out
#SBATCH --error=output/pipeline_%A_%a.err
#SBATCH --time=6:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=60
#SBATCH --partition=cpu

# Prokka output directory (updated for new run with 8,596 genomes)
PROKKA_DIR="../prokka_output"

# Unified BLAST pipeline runner
# This script can run in different modes:
#   1. Array mode (default): Process genomes in batches for BLAST search
#   2. Single mode: Process results and generate statistics
#
# Usage:
#   sbatch --array=1-86 run_pipeline.sh       # Run array job for BLAST search
#   sbatch run_pipeline.sh process            # Process results only
#   sbatch run_pipeline.sh stats              # Generate statistics only
#   sbatch run_pipeline.sh all                # Run complete pipeline (process + overview + stats)

set -euo pipefail

# Parse command line arguments
MODE="search"  # Default to search mode for array jobs
GENOME_LIST_FILE=""

# First positional argument defines the mode (if provided and not an option)
if [[ $# -gt 0 && "$1" != --* ]]; then
    MODE="$1"
    shift
fi

while [[ $# -gt 0 ]]; do
    case $1 in
        --genome-list)
            if [[ -z "${2:-}" ]]; then
                echo "ERROR: --genome-list requires a file path" >&2
                exit 1
            fi
            GENOME_LIST_FILE="$2"
            shift 2
            ;;
        --help)
            echo "Usage: sbatch [--array=...] run_pipeline.sh [MODE] [--genome-list FILE]"
            echo "Modes: search (default), process, overview, stats, all"
            exit 0
            ;;
        *)
            echo "ERROR: Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

if [[ -n "$GENOME_LIST_FILE" ]]; then
    if [[ ! "$GENOME_LIST_FILE" = /* ]]; then
        if [[ -n "$SLURM_SUBMIT_DIR" ]]; then
            GENOME_LIST_FILE="$SLURM_SUBMIT_DIR/$GENOME_LIST_FILE"
        else
            GENOME_LIST_FILE="$(pwd)/$GENOME_LIST_FILE"
        fi
    fi
    if [[ ! -f "$GENOME_LIST_FILE" ]]; then
        echo "ERROR: genome list file not found: $GENOME_LIST_FILE" >&2
        exit 1
    fi
fi

# Initialize conda
eval "$(/home/pmuench/miniconda3/bin/conda shell.bash hook)"
conda activate efs_diversity
export PATH="/home/pmuench/miniconda3/envs/efs_diversity/bin:$PATH"
hash -r

# Set Python path to ensure imports work
export PYTHONPATH="${PYTHONPATH:-}:$(pwd)"

echo "=========================================="
echo "BLAST Pipeline Execution"
echo "Mode: $MODE"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID:-N/A}"
echo "Date: $(date)"
echo "Working directory: $(pwd)"
echo "=========================================="

# Function to run the pipeline
run_pipeline() {
    local mode=$1
    local batch_id=${2:-}
    
    if [ -n "$batch_id" ]; then
        # Running as array job for BLAST search
        echo "Running pipeline in $mode mode for batch $batch_id"
        cmd=(python blast_pipeline.py \
            --mode "$mode" \
            --batch-id "$batch_id" \
            --batch-size 100 \
            --threads "$SLURM_CPUS_PER_TASK" \
            --prokka-dir "$PROKKA_DIR" \
            --blast-modes coding_protein coding_nt prokka_variants noncoding)
        if [[ -n "$GENOME_LIST_FILE" ]]; then
            cmd+=(--genome-list "$GENOME_LIST_FILE")
        fi
        "${cmd[@]}"
    else
        # Running in single mode for processing/stats
        echo "Running pipeline in $mode mode (single job)"
        cmd=(python blast_pipeline.py \
            --mode "$mode" \
            --threads "$SLURM_CPUS_PER_TASK" \
            --prokka-dir "$PROKKA_DIR")
        if [[ -n "$GENOME_LIST_FILE" ]]; then
            cmd+=(--genome-list "$GENOME_LIST_FILE")
        fi
        "${cmd[@]}"
    fi
}

# Main execution logic
case "$MODE" in
    search)
        # BLAST search mode - requires array job
        if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
            echo "ERROR: search mode requires array job submission"
            echo "Use: sbatch --array=1-86 run_pipeline.sh"
            exit 1
        fi
        run_pipeline "search" "$SLURM_ARRAY_TASK_ID"
        ;;
        
    process)
        # Process BLAST results
        run_pipeline "process"
        ;;
        
    overview)
        # Create BLAST overview
        run_pipeline "overview"
        ;;
        
    stats)
        # Generate manuscript statistics
        run_pipeline "stats"
        ;;
        
    all)
        # Run complete pipeline
        if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
            # If running as array, do search for this batch
            run_pipeline "search" "$SLURM_ARRAY_TASK_ID"
        else
            # If single job, do processing and stats
            echo "Running complete analysis pipeline (process + overview + stats)"
            run_pipeline "process"
            run_pipeline "overview"
            run_pipeline "stats"
        fi
        ;;
        
    *)
        echo "ERROR: Unknown mode: $MODE"
        echo "Valid modes: search, process, overview, stats, all"
        exit 1
        ;;
esac

echo "=========================================="
echo "Pipeline execution completed at $(date)"
echo "=========================================="

# If this was the last array job, suggest next steps
if [ "$MODE" = "search" ] && [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    if [ "$SLURM_ARRAY_TASK_ID" = "86" ]; then
        echo ""
        echo "This appears to be the last batch."
        echo "After all array jobs complete, run:"
        echo "  sbatch run_pipeline.sh all"
        echo ""
        echo "Or run individual steps:"
        echo "  sbatch run_pipeline.sh process"
        echo "  sbatch run_pipeline.sh overview"
        echo "  sbatch run_pipeline.sh stats"
        echo ""
    fi
fi
