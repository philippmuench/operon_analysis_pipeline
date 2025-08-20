#!/bin/bash
#SBATCH --job-name=blast_manuscript_stats
#SBATCH --output=manuscript_stats_%j.out
#SBATCH --error=manuscript_stats_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --partition=cpu

echo "=========================================="
echo "Starting BLAST manuscript statistics generation"
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

# Check reference sequences
echo "Checking for reference query sequences..."
REF_DIR="../02_reference_operon_extraction/output"
if [ -d "$REF_DIR" ]; then
    echo "Found reference directory: $REF_DIR"
    for file in operon_genes_protein.fasta operon_genes_nt.fasta operon_noncoding_nt.fasta; do
        if [ -f "$REF_DIR/$file" ]; then
            SEQ_COUNT=$(grep -c "^>" "$REF_DIR/$file" 2>/dev/null || echo "0")
            echo "  $file: $SEQ_COUNT sequences"
        fi
    done
else
    echo "WARNING: Reference directory not found at $REF_DIR"
fi

# Check BLAST output directories
echo ""
echo "Checking for BLAST output directories..."
for dir in output/blast_results output/blast_results_nt output/blast_results_prokka_variants; do
    if [ -d "$dir" ]; then
        FILE_COUNT=$(ls -1 "$dir"/*.txt 2>/dev/null | wc -l)
        echo "  $dir: $FILE_COUNT result files"
    else
        echo "  $dir: Directory not found"
    fi
done

# Check for summary files
echo ""
echo "Checking for summary files..."
if [ -f "output/operon_simple_summary.csv" ]; then
    LINE_COUNT=$(wc -l < output/operon_simple_summary.csv)
    echo "  operon_simple_summary.csv: $LINE_COUNT lines"
fi

echo ""
echo "=========================================="
echo "Running manuscript_numbers.py..."
echo "This will analyze BLAST results across all search strategies"
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
    echo "Last 15 lines of output (summary):"
    echo "----------------------------------------"
    tail -15 manuscript_stats.txt
fi