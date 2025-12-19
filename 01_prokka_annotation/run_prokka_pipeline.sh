#!/bin/bash
#SBATCH --job-name=prokka_pipeline
#SBATCH --output=output/prokka_%A_%a.out
#SBATCH --error=output/prokka_%A_%a.err
#SBATCH --array=1-86%20
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=cpu

# Unified Prokka Pipeline SLURM Script
# Usage:
#   sbatch run_prokka_pipeline.sh                  # Run full pipeline
#   sbatch run_prokka_pipeline.sh --test          # Run test mode (50 genomes)
#   sbatch run_prokka_pipeline.sh --validate      # Only validate existing outputs
#   sbatch run_prokka_pipeline.sh --rerun         # Rerun failed genomes
#   sbatch run_prokka_pipeline.sh --check         # Check progress
#   sbatch run_prokka_pipeline.sh --stats         # Generate manuscript statistics

# Parse command line arguments
MODE="run"
TEST_MODE=""
GENOME_LIST_FILE=""
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --test) TEST_MODE="--test" ;;
        --validate) MODE="validate" ;;
        --rerun) MODE="rerun" ;;
        --check) MODE="check" ;;
        --stats) MODE="stats" ;;
        --genome-list)
            if [[ -z "$2" ]]; then
                echo "Error: --genome-list requires a file path"
                exit 1
            fi
            GENOME_LIST_FILE="$2"
            shift 2
            continue
            ;;
        --help) 
            echo "Usage: sbatch run_prokka_pipeline.sh [OPTIONS]"
            echo "Options:"
            echo "  --test      Run in test mode (first 50 genomes)"
            echo "  --validate  Only validate existing outputs"
            echo "  --rerun     Rerun failed genomes"
            echo "  --check     Check progress"
            echo "  --stats     Generate manuscript statistics"
            echo "  --genome-list FILE  Use only the genomes listed in FILE"
            exit 0
            ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

if [[ -n "$GENOME_LIST_FILE" ]]; then
    if [[ ! "$GENOME_LIST_FILE" = /* ]]; then
        if [[ -n "$SLURM_SUBMIT_DIR" ]]; then
            GENOME_LIST_FILE="$SLURM_SUBMIT_DIR/$GENOME_LIST_FILE"
        else
            GENOME_LIST_FILE="$(pwd)/$GENOME_LIST_FILE"
        fi
    fi
fi

# Initialize conda
eval "$(/home/pmuench/miniconda3/bin/conda shell.bash hook)"
conda activate efs_diversity
export PATH="/home/pmuench/miniconda3/envs/efs_diversity/bin:$PATH"
hash -r

# Ensure output directory exists
mkdir -p output

# Get script directory (handle SLURM copying the script to a spool directory)
SCRIPT_DIR=""
if [[ -n "$SLURM_SUBMIT_DIR" ]]; then
    for candidate in "$SLURM_SUBMIT_DIR" "$SLURM_SUBMIT_DIR/01_prokka_annotation"; do
        if [[ -f "$candidate/prokka_pipeline.py" ]]; then
            SCRIPT_DIR="$candidate"
            break
        fi
    done
fi

if [[ -z "$SCRIPT_DIR" ]]; then
    SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
fi

if [[ ! -f "$SCRIPT_DIR/prokka_pipeline.py" ]]; then
    echo "Error: could not locate prokka_pipeline.py from $SCRIPT_DIR."
    echo "Please submit the job from the 01_prokka_annotation directory or set --chdir."
    exit 1
fi

cd "$SCRIPT_DIR"

if [[ -n "$GENOME_LIST_FILE" ]] && [[ ! -f "$GENOME_LIST_FILE" ]]; then
    echo "Error: genome list file not found: $GENOME_LIST_FILE"
    exit 1
fi

PYTHON_EXTRA_ARGS=()
if [[ -n "$TEST_MODE" ]]; then
    PYTHON_EXTRA_ARGS+=("$TEST_MODE")
fi
if [[ -n "$GENOME_LIST_FILE" ]]; then
    PYTHON_EXTRA_ARGS+=(--genome-list "$GENOME_LIST_FILE")
fi

# Make Python script executable
chmod +x prokka_pipeline.py

# Execute based on mode
case $MODE in
    run)
        # For run mode, handle array job
        if [ -n "$TEST_MODE" ] && [ "$SLURM_ARRAY_TASK_ID" -gt 1 ]; then
            echo "Test mode only needs array task 1. Exiting array task $SLURM_ARRAY_TASK_ID"
            exit 0
        fi
        
        # First prepare genome list if it doesn't exist
        if [ -z "$GENOME_LIST_FILE" ] && [ ! -f "genome_list.txt" ] && [ -z "$TEST_MODE" ]; then
            echo "Preparing genome list..."
            python prokka_pipeline.py prepare
        fi
        if [ -z "$GENOME_LIST_FILE" ] && [ -n "$TEST_MODE" ] && [ ! -f "genome_list_test.txt" ]; then
            echo "Preparing test genome list..."
            python prokka_pipeline.py prepare $TEST_MODE
        fi
        
        # Run Prokka batch for this array task
        echo "Running Prokka for array task $SLURM_ARRAY_TASK_ID"
        CMD=(python prokka_pipeline.py run \
            --array-task-id "$SLURM_ARRAY_TASK_ID" \
            --array-task-max "$SLURM_ARRAY_TASK_MAX")
        if [[ ${#PYTHON_EXTRA_ARGS[@]} -gt 0 ]]; then
            CMD+=("${PYTHON_EXTRA_ARGS[@]}")
        fi
        "${CMD[@]}"
        ;;
    
    validate)
        # Validation doesn't need array job
        if [ "$SLURM_ARRAY_TASK_ID" -ne 1 ]; then
            echo "Validation only needs to run once. Exiting array task $SLURM_ARRAY_TASK_ID"
            exit 0
        fi
        echo "Validating Prokka outputs..."
        python prokka_pipeline.py validate "${PYTHON_EXTRA_ARGS[@]}"
        ;;
    
    rerun)
        # Rerun doesn't use array job
        if [ "$SLURM_ARRAY_TASK_ID" -ne 1 ]; then
            echo "Rerun only needs to run once. Exiting array task $SLURM_ARRAY_TASK_ID"
            exit 0
        fi
        echo "Rerunning failed genomes..."
        python prokka_pipeline.py rerun "${PYTHON_EXTRA_ARGS[@]}"
        
        # After rerun, validate again
        echo "Validating after rerun..."
        python prokka_pipeline.py validate "${PYTHON_EXTRA_ARGS[@]}"
        ;;
    
    check)
        # Check doesn't need array job
        if [ "$SLURM_ARRAY_TASK_ID" -ne 1 ]; then
            echo "Check only needs to run once. Exiting array task $SLURM_ARRAY_TASK_ID"
            exit 0
        fi
        echo "Checking Prokka progress..."
        python prokka_pipeline.py check "${PYTHON_EXTRA_ARGS[@]}"
        ;;
    
    stats)
        # Stats doesn't need array job
        if [ "$SLURM_ARRAY_TASK_ID" -ne 1 ]; then
            echo "Stats only needs to run once. Exiting array task $SLURM_ARRAY_TASK_ID"
            exit 0
        fi
        echo "Generating manuscript statistics..."
        python prokka_pipeline.py stats "${PYTHON_EXTRA_ARGS[@]}"
        ;;
esac

echo "Job completed at: $(date)"
