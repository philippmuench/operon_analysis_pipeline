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
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --test) TEST_MODE="--test" ;;
        --validate) MODE="validate" ;;
        --rerun) MODE="rerun" ;;
        --check) MODE="check" ;;
        --stats) MODE="stats" ;;
        --help) 
            echo "Usage: sbatch run_prokka_pipeline.sh [OPTIONS]"
            echo "Options:"
            echo "  --test      Run in test mode (first 50 genomes)"
            echo "  --validate  Only validate existing outputs"
            echo "  --rerun     Rerun failed genomes"
            echo "  --check     Check progress"
            echo "  --stats     Generate manuscript statistics"
            exit 0
            ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Initialize conda
eval "$(/home/pmuench/miniconda3/bin/conda shell.bash hook)"
conda activate efs_diversity

# Ensure output directory exists
mkdir -p output

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

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
        if [ ! -f "genome_list.txt" ] && [ -z "$TEST_MODE" ]; then
            echo "Preparing genome list..."
            python prokka_pipeline.py prepare
        fi
        if [ -n "$TEST_MODE" ] && [ ! -f "genome_list_test.txt" ]; then
            echo "Preparing test genome list..."
            python prokka_pipeline.py prepare $TEST_MODE
        fi
        
        # Run Prokka batch for this array task
        echo "Running Prokka for array task $SLURM_ARRAY_TASK_ID"
        python prokka_pipeline.py run \
            --array-task-id $SLURM_ARRAY_TASK_ID \
            --array-task-max $SLURM_ARRAY_TASK_MAX \
            $TEST_MODE
        ;;
    
    validate)
        # Validation doesn't need array job
        if [ "$SLURM_ARRAY_TASK_ID" -ne 1 ]; then
            echo "Validation only needs to run once. Exiting array task $SLURM_ARRAY_TASK_ID"
            exit 0
        fi
        echo "Validating Prokka outputs..."
        python prokka_pipeline.py validate $TEST_MODE
        ;;
    
    rerun)
        # Rerun doesn't use array job
        if [ "$SLURM_ARRAY_TASK_ID" -ne 1 ]; then
            echo "Rerun only needs to run once. Exiting array task $SLURM_ARRAY_TASK_ID"
            exit 0
        fi
        echo "Rerunning failed genomes..."
        python prokka_pipeline.py rerun $TEST_MODE
        
        # After rerun, validate again
        echo "Validating after rerun..."
        python prokka_pipeline.py validate $TEST_MODE
        ;;
    
    check)
        # Check doesn't need array job
        if [ "$SLURM_ARRAY_TASK_ID" -ne 1 ]; then
            echo "Check only needs to run once. Exiting array task $SLURM_ARRAY_TASK_ID"
            exit 0
        fi
        echo "Checking Prokka progress..."
        python prokka_pipeline.py check $TEST_MODE
        ;;
    
    stats)
        # Stats doesn't need array job
        if [ "$SLURM_ARRAY_TASK_ID" -ne 1 ]; then
            echo "Stats only needs to run once. Exiting array task $SLURM_ARRAY_TASK_ID"
            exit 0
        fi
        echo "Generating manuscript statistics..."
        python prokka_pipeline.py stats $TEST_MODE
        ;;
esac

echo "Job completed at: $(date)"