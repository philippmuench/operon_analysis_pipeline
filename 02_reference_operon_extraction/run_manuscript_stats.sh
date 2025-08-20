#!/bin/bash
#SBATCH --job-name=ref_operon_stats
#SBATCH --output=manuscript_stats_%j.out
#SBATCH --error=manuscript_stats_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu

echo "=========================================="
echo "Starting manuscript statistics generation"
echo "Job ID: $SLURM_JOB_ID"
echo "Date: $(date)"
echo "=========================================="

# Initialize conda
echo "Activating conda environment..."
eval "$(/home/pmuench/miniconda3/bin/conda shell.bash hook)"
conda activate efs_diversity
echo "Conda environment activated: $CONDA_DEFAULT_ENV"

# Show working directory
echo "Working directory: $(pwd)"

# Check input files
echo "Checking for input files..."
if [ -f "operon.gb" ]; then
    echo "Found GenBank file: operon.gb"
    echo "File size: $(ls -lh operon.gb | awk '{print $5}')"
else
    echo "WARNING: GenBank file 'operon.gb' not found"
fi

# Check output directory
echo "Checking for output files..."
if [ -d "output" ]; then
    FILE_COUNT=$(ls -1 output/*.fasta output/*.tsv 2>/dev/null | wc -l)
    echo "Found $FILE_COUNT output files in output/"
else
    echo "WARNING: Output directory not found"
fi

echo "=========================================="
echo "Running manuscript_numbers.py..."
echo "=========================================="

# Run with verbose Python output
python -u manuscript_numbers.py manuscript_stats.txt

echo "=========================================="
echo "Statistics generation complete at $(date)"
echo "Output saved to manuscript_stats.txt"
echo "=========================================="

# Show summary of output file
if [ -f "manuscript_stats.txt" ]; then
    echo "Output file size: $(ls -lh manuscript_stats.txt | awk '{print $5}')"
    echo ""
    echo "First 20 lines of output:"
    head -20 manuscript_stats.txt
fi