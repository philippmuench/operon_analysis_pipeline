#!/bin/bash
#SBATCH --job-name=dnds_pipeline
#SBATCH --output=dnds_pipeline_%j.out
#SBATCH --error=dnds_pipeline_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=20
#SBATCH --partition=cpu

echo "=========================================="
echo "Starting Consolidated dN/dS Analysis Pipeline"
echo "Job ID: $SLURM_JOB_ID"
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "Working directory: $(pwd)"
echo "=========================================="

# Initialize conda
echo ""
echo "Setting up environment..."
eval "$(/home/pmuench/miniconda3/bin/conda shell.bash hook)"
conda activate efs_diversity
echo "Conda environment activated: $CONDA_DEFAULT_ENV"

# Check Python version and required packages
echo ""
echo "Environment check:"
python --version
echo "Checking required packages..."
python -c "import Bio; print(f'  BioPython: {Bio.__version__}')" 2>/dev/null || echo "  BioPython: Not found"
python -c "import pandas; print(f'  Pandas: {pandas.__version__}')" 2>/dev/null || echo "  Pandas: Not found"
python -c "import numpy; print(f'  NumPy: {numpy.__version__}')" 2>/dev/null || echo "  NumPy: Not found"
python -c "import matplotlib; print(f'  Matplotlib: {matplotlib.__version__}')" 2>/dev/null || echo "  Matplotlib: Not found"
python -c "import seaborn; print(f'  Seaborn: {seaborn.__version__}')" 2>/dev/null || echo "  Seaborn: Not found"

# Set input directories
OPERON_DIR="../05_operon_assembly_extraction/output/msa/dna_alignments"
CORE_DIR="../04_core_gene_analysis/output/core_gene_alignments"
OUTPUT_DIR="output"

# Check for input files
echo ""
echo "=========================================="
echo "Checking input directories..."
echo "=========================================="

# Check operon alignments
if [ -d "$OPERON_DIR" ]; then
    OPERON_COUNT=$(find "$OPERON_DIR" -type f \( -name "*.fasta" -o -name "*_aligned.fasta" \) 2>/dev/null | wc -l)
    echo "✓ Found operon alignment directory: $OPERON_DIR"
    echo "  Contains $OPERON_COUNT alignment files"
    if [ "$OPERON_COUNT" -eq 0 ]; then
        echo "  ⚠ WARNING: No alignment files found in operon directory"
        echo "  Please run step 05 first: cd ../05_operon_assembly_extraction && sbatch run_operon_extraction.sh"
    fi
else
    echo "✗ Operon alignment directory not found: $OPERON_DIR"
    echo "  Alternative paths to check:"
    echo "    - ../05_operon_assembly_extraction/output/mappings/aa_nt_mapping/prokka/msa/"
    echo "    - ../05_operon_assembly_extraction/output/mappings/nt_nt_mapping/prokka_genome/msa/"
    OPERON_DIR=""
fi

# Check core gene alignments
if [ -d "$CORE_DIR" ]; then
    CORE_COUNT=$(find "$CORE_DIR" -type f \( -name "*.fasta" -o -name "*_aligned.fasta" \) 2>/dev/null | wc -l)
    echo "✓ Found core gene alignment directory: $CORE_DIR"
    echo "  Contains $CORE_COUNT alignment files"
    if [ "$CORE_COUNT" -eq 0 ]; then
        echo "  ⚠ WARNING: No alignment files found in core gene directory"
        echo "  Please run step 04 first: cd ../04_core_gene_analysis && sbatch run_core_analysis.sh"
    fi
else
    echo "✗ Core gene alignment directory not found: $CORE_DIR"
    echo "  Please run step 04 first to generate core gene alignments"
    CORE_DIR=""
fi

# Check if we have at least one valid input directory
if [ -z "$OPERON_DIR" ] && [ -z "$CORE_DIR" ]; then
    echo ""
    echo "ERROR: No valid input directories found!"
    echo "Please ensure previous pipeline steps have been completed."
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/plots"
mkdir -p "$OUTPUT_DIR/tables"
mkdir -p "$OUTPUT_DIR/logs"

echo ""
echo "=========================================="
echo "Running dN/dS Analysis Pipeline..."
echo "=========================================="
echo "Configuration:"
echo "  Operon alignments: ${OPERON_DIR:-'Not available'}"
echo "  Core alignments: ${CORE_DIR:-'Not available'}"
echo "  Output directory: $OUTPUT_DIR"
echo "  Threads: $SLURM_CPUS_PER_TASK"
echo ""

# Build command based on available inputs
CMD="python dnds_pipeline.py --output-dir $OUTPUT_DIR --threads $SLURM_CPUS_PER_TASK"

if [ -n "$OPERON_DIR" ]; then
    CMD="$CMD --operon-dir $OPERON_DIR"
fi

if [ -n "$CORE_DIR" ]; then
    CMD="$CMD --core-dir $CORE_DIR"
fi

# Run the pipeline
echo "Executing: $CMD"
echo ""
$CMD

# Check if pipeline completed successfully
if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Pipeline completed successfully!"
    echo "=========================================="
    
    # Show output summary
    echo ""
    echo "Output Summary:"
    echo "---------------"
    
    # List generated files
    if [ -d "$OUTPUT_DIR" ]; then
        echo ""
        echo "Generated files:"
        
        # Tables
        if [ -d "$OUTPUT_DIR/tables" ] && [ "$(ls -A $OUTPUT_DIR/tables 2>/dev/null)" ]; then
            echo "  Tables:"
            ls -lh "$OUTPUT_DIR/tables"/*.csv 2>/dev/null | awk '{print "    - " $9 ": " $5}'
        fi
        
        # Plots
        if [ -d "$OUTPUT_DIR/plots" ] && [ "$(ls -A $OUTPUT_DIR/plots 2>/dev/null)" ]; then
            echo "  Plots:"
            ls -lh "$OUTPUT_DIR/plots"/*.{png,pdf} 2>/dev/null | awk '{print "    - " $9 ": " $5}'
        fi
        
        # Reports
        if [ -f "$OUTPUT_DIR/dnds_summary_report.txt" ]; then
            echo "  Reports:"
            ls -lh "$OUTPUT_DIR"/*.txt 2>/dev/null | awk '{print "    - " $9 ": " $5}'
        fi
        
        # Show quick summary if available
        if [ -f "$OUTPUT_DIR/dnds_summary_statistics.txt" ]; then
            echo ""
            echo "Quick Summary Statistics:"
            echo "-------------------------"
            cat "$OUTPUT_DIR/dnds_summary_statistics.txt"
            echo "-------------------------"
        fi
        
        # Show first few lines of the summary report
        if [ -f "$OUTPUT_DIR/dnds_summary_report.txt" ]; then
            echo ""
            echo "Report Preview (first 30 lines):"
            echo "--------------------------------"
            head -30 "$OUTPUT_DIR/dnds_summary_report.txt"
            echo "..."
            echo "--------------------------------"
            echo "Full report: $OUTPUT_DIR/dnds_summary_report.txt"
        fi
    fi
    
else
    echo ""
    echo "=========================================="
    echo "ERROR: Pipeline failed!"
    echo "=========================================="
    echo "Please check the error messages above and the log files."
    exit 1
fi

echo ""
echo "=========================================="
echo "Job completed at $(date)"
echo "Total runtime: $SECONDS seconds"
echo "=========================================="