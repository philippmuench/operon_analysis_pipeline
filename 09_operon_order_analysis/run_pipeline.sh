#!/bin/bash
#SBATCH --job-name=operon_order_pipeline
#SBATCH --output=pipeline_%j.out
#SBATCH --error=pipeline_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu

echo "=========================================="
echo "Operon Order Analysis Pipeline"
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
echo ""

# Check for input BLAST results
echo "Checking for input files..."
BLAST_DIR="../03_blast_search/output/blast_results"

if [ -d "$BLAST_DIR" ]; then
    FILE_COUNT=$(find "$BLAST_DIR" -name "ENT_*_genes_blast.txt" | wc -l)
    echo "Found $FILE_COUNT genome BLAST result files"
else
    echo "ERROR: BLAST results directory not found at $BLAST_DIR"
    echo "Please run step 03 first: cd ../03_blast_search && sbatch run_blast_search.sh"
    exit 1
fi

echo ""
echo "=========================================="
echo "Running Operon Order Analysis Pipeline"
echo "=========================================="

# Set metadata file path
METADATA_FILE="../00_annotation/8587_Efs_metadata_ASbarcode.txt"

echo "Will also perform metadata stratification using $METADATA_FILE"
echo ""

# Run the unified pipeline with all steps
python operon_order_pipeline.py \
    --blast-dir "$BLAST_DIR" \
    --output-dir output \
    --min-identity 80 \
    --metadata "$METADATA_FILE" \
    --steps all

# Check if pipeline succeeded
if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Pipeline completed successfully!"
    echo "=========================================="
    
    # Show output summary
    OUTPUT_DIR="output"
    if [ -d "$OUTPUT_DIR" ]; then
        echo ""
        echo "Generated files:"
        ls -lh "$OUTPUT_DIR"/*.csv "$OUTPUT_DIR"/*.txt "$OUTPUT_DIR"/*.png "$OUTPUT_DIR"/*.pdf 2>/dev/null | awk '{print "  " $9 ": " $5}'
        
        # Display manuscript statistics if available
        if [ -f "$OUTPUT_DIR/manuscript_stats.txt" ]; then
            echo ""
            echo "=========================================="
            echo "Manuscript Statistics:"
            echo "=========================================="
            head -30 "$OUTPUT_DIR/manuscript_stats.txt"
            echo "..."
            echo "(Full statistics in $OUTPUT_DIR/manuscript_stats.txt)"
        fi
        
        # Display summary statistics
        if [ -f "$OUTPUT_DIR/operon_order_summary.txt" ]; then
            echo ""
            echo "=========================================="
            echo "Analysis Summary:"
            echo "=========================================="
            head -20 "$OUTPUT_DIR/operon_order_summary.txt"
        fi
        
        # Display metadata stratification if available
        if [ -f "$OUTPUT_DIR/operon_order_metadata_stratification.txt" ]; then
            echo ""
            echo "=========================================="
            echo "Metadata Stratification Summary:"
            echo "=========================================="
            cat "$OUTPUT_DIR/operon_order_metadata_stratification.txt"
            echo ""
            echo "Additional stratification files generated:"
            echo "  - operon_order_by_source_niche.tsv"
            echo "  - operon_order_by_country.tsv"
            echo "  - gene_presence_by_source_niche.tsv"
        fi
    fi
else
    echo ""
    echo "ERROR: Pipeline failed!"
    echo "Check the error messages above for details."
    exit 1
fi

echo ""
echo "=========================================="
echo "Analysis complete at $(date)"
echo "=========================================="