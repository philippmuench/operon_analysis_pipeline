#!/bin/bash
#SBATCH --job-name=manuscript_stats
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

# Check Prokka version
echo "Checking Prokka version..."
prokka --version 2>&1 || echo "Could not get Prokka version directly"

# Show working directory
echo "Working directory: $(pwd)"

# Check input directory
echo "Checking for input genomes..."
GENOME_DIR="../../Efs_assemblies"
if [ -d "$GENOME_DIR" ]; then
    GENOME_COUNT=$(ls -1 "$GENOME_DIR"/*.fasta.gz 2>/dev/null | wc -l)
    echo "Found $GENOME_COUNT genome files in $GENOME_DIR"
else
    echo "WARNING: Genome directory not found at $GENOME_DIR"
fi

# Check Prokka output directory
echo "Checking for Prokka outputs..."
PROKKA_DIR="../prokka_output"
if [ -d "$PROKKA_DIR" ]; then
    DIR_COUNT=$(find "$PROKKA_DIR" -maxdepth 1 -type d | wc -l)
    echo "Found $((DIR_COUNT-1)) directories in $PROKKA_DIR"
else
    PROKKA_DIR="output/prokka_results"
    if [ -d "$PROKKA_DIR" ]; then
        DIR_COUNT=$(find "$PROKKA_DIR" -maxdepth 1 -type d | wc -l)
        echo "Using test directory: Found $((DIR_COUNT-1)) directories in $PROKKA_DIR"
    else
        echo "WARNING: No Prokka output directory found"
    fi
fi

echo "=========================================="
echo "Running manuscript_numbers.py..."
echo "This may take 10-15 minutes for 8,587 genomes"
echo "The script is analyzing:"
echo "  - Counting all GFF files for gene statistics"
echo "  - Counting all FAA files for protein statistics"
echo "  - Calculating averages and outliers"
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
    echo "First 10 lines of output:"
    head -10 manuscript_stats.txt
fi