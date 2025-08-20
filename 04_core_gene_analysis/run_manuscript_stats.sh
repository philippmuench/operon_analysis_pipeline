#!/bin/bash
#SBATCH --job-name=core_gene_manuscript_stats
#SBATCH --output=manuscript_stats_%j.out
#SBATCH --error=manuscript_stats_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu

echo "=========================================="
echo "Starting Core Gene Analysis manuscript statistics generation"
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

# Check for input data
echo ""
echo "Checking for input data..."

# Check Prokka results
PROKKA_DIR="../01_prokka_annotation/output/prokka_results"
if [ -d "$PROKKA_DIR" ]; then
    GENOME_COUNT=$(find "$PROKKA_DIR" -maxdepth 1 -type d | wc -l)
    echo "  Prokka results: $((GENOME_COUNT-1)) genome directories"
else
    echo "  WARNING: Prokka results directory not found at $PROKKA_DIR"
fi

# Check for core gene analysis outputs
echo ""
echo "Checking for analysis outputs..."

# Check core genes list
if [ -f "output/core_genes_95pct.txt" ]; then
    CORE_COUNT=$(wc -l < output/core_genes_95pct.txt)
    echo "  Core genes list: $CORE_COUNT genes"
else
    echo "  Core genes list: Not found"
fi

# Check prevalence statistics
if [ -f "output/gene_prevalence_stats.csv" ]; then
    LINE_COUNT=$(wc -l < output/gene_prevalence_stats.csv)
    echo "  Gene prevalence stats: $((LINE_COUNT-1)) genes"
else
    echo "  Gene prevalence stats: Not found"
fi

# Check sequence directories
if [ -d "output/core_gene_sequences" ]; then
    SEQ_COUNT=$(ls -1 output/core_gene_sequences/*.fasta 2>/dev/null | wc -l)
    echo "  Core gene sequences: $SEQ_COUNT FASTA files"
fi

if [ -d "output/core_gene_alignments" ]; then
    ALIGN_COUNT=$(ls -1 output/core_gene_alignments/*_aligned.fasta 2>/dev/null | wc -l)
    echo "  Core gene alignments: $ALIGN_COUNT alignment files"
fi

# Check conservation metrics
if [ -f "output/core_gene_conservation_metrics.csv" ]; then
    CONS_COUNT=$(wc -l < output/core_gene_conservation_metrics.csv)
    echo "  Conservation metrics: $((CONS_COUNT-1)) genes analyzed"
fi

echo ""
echo "=========================================="
echo "Running manuscript_numbers.py..."
echo "This will analyze core gene statistics"
echo "=========================================="

# Run with verbose Python output and save to file
python -u manuscript_numbers.py manuscript_stats.txt

echo ""
echo "=========================================="
echo "Statistics generation complete at $(date)"
echo "Output saved to manuscript_stats.txt"
echo "=========================================="

# Show summary of output file
if [ -f "manuscript_stats.txt" ]; then
    echo ""
    echo "Output file size: $(ls -lh manuscript_stats.txt | awk '{print $5}')"
    echo ""
    echo "First 20 lines of output:"
    echo "----------------------------------------"
    head -20 manuscript_stats.txt
    echo "----------------------------------------"
    echo ""
    echo "Last 20 lines of output (summary):"
    echo "----------------------------------------"
    tail -20 manuscript_stats.txt
fi