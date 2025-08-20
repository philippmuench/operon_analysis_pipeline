#!/bin/bash
#SBATCH --job-name=operon_order_manuscript_stats
#SBATCH --output=manuscript_stats_%j.out
#SBATCH --error=manuscript_stats_%j.err
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu

echo "=========================================="
echo "Generating Manuscript Statistics"
echo "Job ID: $SLURM_JOB_ID"
echo "Date: $(date)"
echo "=========================================="

# Initialize conda
echo "Activating conda environment..."
eval "$(/home/pmuench/miniconda3/bin/conda shell.bash hook)"
conda activate efs_diversity
echo "Conda environment activated: $CONDA_DEFAULT_ENV"

# Check for results file
if [ ! -f "output/operon_order_analysis.csv" ]; then
    echo "ERROR: Analysis results not found at output/operon_order_analysis.csv"
    echo "Please run the operon order analysis first:"
    echo "  sbatch run_operon_order_analysis.sh"
    exit 1
fi

# Generate manuscript statistics
echo ""
echo "Generating statistics from operon order analysis..."
python manuscript_numbers.py --output output/manuscript_stats.txt

# Display the results
if [ -f "output/manuscript_stats.txt" ]; then
    echo ""
    echo "=========================================="
    echo "Manuscript Statistics:"
    echo "=========================================="
    cat output/manuscript_stats.txt
    echo ""
    echo "=========================================="
    echo "Statistics saved to: output/manuscript_stats.txt"
else
    echo "ERROR: Failed to generate manuscript statistics"
    exit 1
fi

echo ""
echo "Complete at $(date)"