#!/bin/bash
#SBATCH --job-name=core_gene_pipeline
#SBATCH --output=pipeline_%j.out
#SBATCH --error=pipeline_%j.err
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=90
#SBATCH --mem=300G
#SBATCH --partition=cpu

# Core Gene Analysis Pipeline SLURM Script
# =========================================
# This script runs the consolidated core gene analysis pipeline

# Parse command line arguments
START_STEP=1
THREADS=${SLURM_CPUS_PER_TASK:-90}
PROKKA_DIR="../01_prokka_annotation/output/prokka_results"
OUTPUT_DIR="output"
THRESHOLD=0.95
MAFFT_TIMEOUT=0
HELP=false

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Core Gene Analysis Pipeline - SLURM submission script"
    echo ""
    echo "Options:"
    echo "  --start-step STEP       Start pipeline from specific step (1-5, default: 1)"
    echo "  --prokka-dir DIR        Directory with Prokka outputs (default: ../01_prokka_annotation/output/prokka_results)"
    echo "  --output-dir DIR        Output directory (default: output)"
    echo "  --threshold FLOAT       Core gene prevalence threshold (default: 0.95)"
    echo "  --threads NUM           Number of threads (default: SLURM allocation or 90)"
    echo "  --mafft-timeout SEC     MAFFT timeout per gene in seconds (default: 0=none)"
    echo "  --help                  Show this help message"
    echo ""
    echo "Pipeline steps:"
    echo "  1. Identify core genes (≥threshold prevalence)"
    echo "  2. Extract sequences from Prokka output"
    echo "  3. Create MSAs using MAFFT"
    echo "  4. Calculate diversity metrics"
    echo "  5. Threshold analysis and plots"
    echo ""
    echo "Examples:"
    echo "  sbatch $0                           # Run complete pipeline"
    echo "  sbatch $0 --start-step 3            # Start from MSA creation"
    echo "  sbatch $0 --threshold 0.99          # Use 99% threshold"
    echo "  sbatch $0 --start-step 5            # Only run threshold analysis"
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --start-step)
            START_STEP="$2"
            shift 2
            ;;
        --prokka-dir)
            PROKKA_DIR="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --threshold)
            THRESHOLD="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --mafft-timeout)
            MAFFT_TIMEOUT="$2"
            shift 2
            ;;
        --help)
            HELP=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

if [ "$HELP" = true ]; then
    usage
    exit 0
fi

# Validate inputs
if ! [[ "$START_STEP" =~ ^[1-5]$ ]]; then
    echo "❌ Error: Invalid start step '$START_STEP'. Must be 1-5."
    usage
    exit 1
fi

# Validate threshold is between 0 and 1
if (( $(echo "$THRESHOLD < 0" | bc -l) )) || (( $(echo "$THRESHOLD > 1" | bc -l) )); then
    echo "❌ Error: Threshold must be between 0 and 1 (got: $THRESHOLD)"
    exit 1
fi

# Set MAFFT temporary directory
export MAFFT_TMPDIR="/vol/tmp"

# Validate MAFFT_TMPDIR
echo "Validating MAFFT configuration..."
if [ -z "$MAFFT_TMPDIR" ]; then
    echo "❌ Error: MAFFT_TMPDIR environment variable is not set"
    exit 1
fi

if [ ! -d "$MAFFT_TMPDIR" ]; then
    echo "❌ Error: MAFFT_TMPDIR directory does not exist: $MAFFT_TMPDIR"
    exit 1
fi

if [ ! -w "$MAFFT_TMPDIR" ]; then
    echo "❌ Error: MAFFT_TMPDIR is not writable: $MAFFT_TMPDIR"
    exit 1
fi

echo "✅ MAFFT_TMPDIR validated: $MAFFT_TMPDIR"

# Print configuration
echo ""
echo "=================================================="
echo "Core Gene Analysis Pipeline"
echo "=================================================="
echo "Started: $(date)"
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Configuration:"
echo "  Prokka directory: $PROKKA_DIR"
echo "  Output directory: $OUTPUT_DIR"
echo "  Threshold: ${THRESHOLD} ($(echo "$THRESHOLD * 100" | bc)%)"
echo "  Threads: $THREADS"
echo "  MAFFT timeout: $MAFFT_TIMEOUT seconds"
echo "  Starting from step: $START_STEP"
echo "=================================================="

# Change to script directory
cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/04_core_gene_analysis

# Run the pipeline
echo ""
echo "Running consolidated pipeline..."
python core_gene_pipeline.py \
    --prokka-dir "$PROKKA_DIR" \
    --output-dir "$OUTPUT_DIR" \
    --threshold "$THRESHOLD" \
    --threads "$THREADS" \
    --start-step "$START_STEP" \
    --mafft-timeout "$MAFFT_TIMEOUT"

PIPELINE_EXIT_CODE=$?

# Check exit status
if [ $PIPELINE_EXIT_CODE -eq 0 ]; then
    echo ""
    echo "=================================================="
    echo "✅ Pipeline completed successfully!"
    echo "Finished: $(date)"
    echo "=================================================="
    
    # Run manuscript statistics if complete pipeline was run
    if [ $START_STEP -eq 1 ] && [ -f "manuscript_numbers.py" ]; then
        echo ""
        echo "Generating manuscript statistics..."
        python manuscript_numbers.py
    fi
    
    exit 0
else
    echo ""
    echo "=================================================="
    echo "❌ Pipeline failed with exit code: $PIPELINE_EXIT_CODE"
    echo "Finished: $(date)"
    echo "=================================================="
    exit $PIPELINE_EXIT_CODE
fi