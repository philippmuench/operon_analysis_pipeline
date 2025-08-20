#!/bin/bash
#SBATCH --job-name=operon_order
#SBATCH --output=operon_order_%j.out
#SBATCH --error=operon_order_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu

echo "=========================================="
echo "Starting Operon Order and Synteny Analysis"
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

# Check for input BLAST results
echo ""
echo "Checking for input files..."

# Check for BLAST results from step 03
BLAST_DIR="../03_blast_search/output/blast_results"
if [ -d "$BLAST_DIR" ]; then
    FILE_COUNT=$(find "$BLAST_DIR" -name "ENT_*_genes_blast.txt" | wc -l)
    echo "Found $FILE_COUNT genome BLAST result files"
    
    # Show a few examples
    echo "Example files:"
    ls "$BLAST_DIR"/ENT_*_genes_blast.txt 2>/dev/null | head -5 | while read file; do
        echo "  $(basename $file)"
    done
else
    echo "ERROR: BLAST results directory not found at $BLAST_DIR"
    echo "   Please run step 03 first: cd ../03_blast_search && sbatch run_blast_search.sh"
    exit 1
fi

echo ""
echo "=========================================="
echo "Analyzing operon gene order and synteny..."
echo "=========================================="

# Run the operon order analysis
python analyze_operon_order.py \
    --blast-dir "$BLAST_DIR" \
    --output-dir output \
    --min-identity 80

echo ""
echo "=========================================="
echo "Creating visualizations..."
echo "=========================================="

# Create visualizations if analysis was successful
if [ -f "output/operon_order_analysis.csv" ]; then
    python visualize_operon_synteny.py
else
    echo "WARNING: Analysis results not found, skipping visualization"
fi

echo ""
echo "=========================================="
echo "Generating manuscript statistics..."
echo "=========================================="

# Generate manuscript statistics if analysis was successful
if [ -f "output/operon_order_analysis.csv" ]; then
    python manuscript_numbers.py --output output/manuscript_stats.txt
    
    # Display the statistics
    if [ -f "output/manuscript_stats.txt" ]; then
        echo ""
        echo "Manuscript Statistics Generated:"
        echo "----------------------------------------"
        cat output/manuscript_stats.txt
        echo "----------------------------------------"
    fi
else
    echo "WARNING: Analysis results not found, skipping manuscript statistics"
fi

echo ""
echo "=========================================="
echo "Analysis complete at $(date)"
echo "=========================================="

# Show output summary
OUTPUT_DIR="output"
if [ -d "$OUTPUT_DIR" ]; then
    echo ""
    echo "Generated files:"
    ls -lh "$OUTPUT_DIR"/*.csv "$OUTPUT_DIR"/*.txt "$OUTPUT_DIR"/*.png "$OUTPUT_DIR"/*.pdf 2>/dev/null | awk '{print "  " $9 ": " $5}'
    
    # Show summary statistics
    if [ -f "$OUTPUT_DIR/operon_order_summary.txt" ]; then
        echo ""
        echo "Operon Order Summary:"
        echo "----------------------------------------"
        head -20 "$OUTPUT_DIR/operon_order_summary.txt"
        echo "----------------------------------------"
    fi
fi