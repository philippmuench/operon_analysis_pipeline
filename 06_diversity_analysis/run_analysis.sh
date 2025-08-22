#!/bin/bash
#SBATCH --job-name=comparative_analysis
#SBATCH --output=comparative_analysis_%j.out
#SBATCH --error=comparative_analysis_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu

echo "=========================================="
echo "Starting Streamlined Comparative Analysis"
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

# Check for input files from previous steps
echo ""
echo "Checking for input files..."

# Check core gene metrics from step 04
CORE_CSV="../04_core_gene_analysis/output/core_gene_conservation_metrics.csv"
if [ -f "$CORE_CSV" ]; then
    LINES=$(wc -l < "$CORE_CSV")
    echo "Core gene metrics found: $((LINES-1)) genes"
else
    echo "WARNING: Core gene metrics not found at $CORE_CSV"
    echo "   Please run step 04 first: cd ../04_core_gene_analysis && sbatch run_core_analysis.sh"
fi

# Check operon metrics from step 05
echo ""
echo "Checking for operon metrics from step 05..."

# Use detailed metrics if available (matches core gene format), otherwise use basic
OPERON_CSV_DETAILED="../05_operon_assembly_extraction/output/operon_conservation_metrics_detailed.csv"
OPERON_CSV_BASIC="../05_operon_assembly_extraction/output/operon_conservation_metrics.csv"

if [ -f "$OPERON_CSV_DETAILED" ]; then
    OPERON_CSV="$OPERON_CSV_DETAILED"
    echo "Using detailed operon conservation metrics (matching core gene format)"
    LINES=$(wc -l < "$OPERON_CSV")
    echo "Operon metrics found: $((LINES-1)) genes"
elif [ -f "$OPERON_CSV_BASIC" ]; then
    OPERON_CSV="$OPERON_CSV_BASIC"
    echo "Using basic operon conservation metrics"
    LINES=$(wc -l < "$OPERON_CSV")
    echo "Operon metrics found: $((LINES-1)) genes"
else
    echo "WARNING: Operon metrics not found at either location"
    
    # Also check legacy directory structure
    OPERON_DIR="../05_operon_assembly_extraction/output/mappings"
    if [ -d "$OPERON_DIR" ]; then
        # Count CSV files
        CSV_COUNT=$(find "$OPERON_DIR" -name "operon_conservation_metrics.csv" | wc -l)
        echo "Found $CSV_COUNT operon conservation metric files in legacy location"
        
        # List strategies
        for csv in $(find "$OPERON_DIR" -name "operon_conservation_metrics.csv"); do
            STRATEGY=$(echo "$csv" | sed 's/.*mappings\///' | sed 's/\/msa.*//')
            LINES=$(wc -l < "$csv")
            echo "   - $STRATEGY: $((LINES-1)) genes"
        done
    else
        echo "ERROR: No operon metrics found"
        echo "   Please run step 05 first: cd ../05_operon_assembly_extraction && sbatch run_operon_extraction.sh"
    fi
fi

echo ""
echo "=========================================="
echo "Running comparative analysis..."
echo "Creating conservation ranking visualizations"
echo "=========================================="

# Parse command line arguments
MODE="comparative"
if [ "$1" == "--plots" ]; then
    MODE="plots"
elif [ "$1" == "--both" ]; then
    MODE="both"
fi

# Run the comparative analysis
if [ "$MODE" == "plots" ]; then
    echo "Creating conservation plots from MSA alignments..."
    python comparative_analysis.py --mode plots
elif [ "$MODE" == "both" ]; then
    echo "Running both comparative analysis and conservation plots..."
    python comparative_analysis.py --mode both
else
    echo "Running comparative analysis only..."
    python comparative_analysis.py --mode comparative
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
    ls -lh "$OUTPUT_DIR"/conservation_*.{png,pdf} "$OUTPUT_DIR"/operon_conservation_summary.csv 2>/dev/null | awk '{print "  " $9 ": " $5}'
    
    # Show summary statistics
    if [ -f "$OUTPUT_DIR/operon_conservation_summary.csv" ]; then
        echo ""
        echo "Operon conservation summary:"
        echo "----------------------------------------"
        cat "$OUTPUT_DIR/operon_conservation_summary.csv"
        echo "----------------------------------------"
    fi
fi