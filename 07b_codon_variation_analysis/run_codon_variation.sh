#!/bin/bash
#SBATCH --job-name=codon_variation
#SBATCH --output=codon_variation_%j.out
#SBATCH --error=codon_variation_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu

echo "=========================================="
echo "Codon Variation-Based Selection Analysis"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "Working directory: $(pwd)"
echo "=========================================="

# Activate conda environment
echo ""
echo "Setting up environment..."
eval "$(/home/pmuench/miniconda3/bin/conda shell.bash hook)"
conda activate efs_diversity
echo "Conda environment activated: $CONDA_DEFAULT_ENV"

# Check Python version
echo ""
echo "Environment check:"
python --version

# Check input files exist
OPERON_FILE="../07_dnds_analysis/output_old/tables/operon_dnds_analysis.csv"
CORE_FILE="../07_dnds_analysis/output_old/tables/core_genes_dnds_analysis.csv"

echo ""
echo "Checking input files..."
if [ -f "$OPERON_FILE" ]; then
    echo "✓ Found operon file: $OPERON_FILE"
    echo "  Size: $(du -h $OPERON_FILE | cut -f1)"
    echo "  Lines: $(wc -l < $OPERON_FILE)"
else
    echo "✗ Operon file not found: $OPERON_FILE"
    echo "  Please run 07_dnds_analysis first"
    exit 1
fi

if [ -f "$CORE_FILE" ]; then
    echo "✓ Found core gene file: $CORE_FILE"
    echo "  Size: $(du -h $CORE_FILE | cut -f1)"
    echo "  Lines: $(wc -l < $CORE_FILE)"
else
    echo "✗ Core gene file not found: $CORE_FILE"
    echo "  Please run 07_dnds_analysis first"
    exit 1
fi

# Create output directory
mkdir -p output/{tables,plots,logs}

echo ""
echo "=========================================="
echo "Running Codon Variation Analysis..."
echo "=========================================="

# Run the analysis
python codon_variation_analysis.py \
    --operon-file "$OPERON_FILE" \
    --core-file "$CORE_FILE" \
    --output-dir output

# Check if analysis completed successfully
if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Analysis completed successfully!"
    echo "=========================================="
    
    # Show output summary
    echo ""
    echo "Output Summary:"
    echo "---------------"
    
    # List generated files
    if [ -d "output" ]; then
        echo ""
        echo "Generated files:"
        
        # Tables
        if [ -d "output/tables" ] && [ "$(ls -A output/tables 2>/dev/null)" ]; then
            echo "  Tables:"
            ls -lh output/tables/*.csv 2>/dev/null | awk '{print "    - " $9 ": " $5}'
        fi
        
        # Plots
        if [ -d "output/plots" ] && [ "$(ls -A output/plots 2>/dev/null)" ]; then
            echo "  Plots:"
            ls -lh output/plots/*.png 2>/dev/null | awk '{print "    - " $9 ": " $5}'
        fi
        
        # Report
        if [ -f "output/codon_variation_report.txt" ]; then
            echo "  Report:"
            ls -lh output/*.txt 2>/dev/null | awk '{print "    - " $9 ": " $5}'
            
            # Show summary from report
            echo ""
            echo "Report Preview:"
            echo "---------------"
            head -50 output/codon_variation_report.txt
            echo "..."
            echo "---------------"
            echo "Full report: output/codon_variation_report.txt"
        fi
    fi
else
    echo ""
    echo "=========================================="
    echo "ERROR: Analysis failed!"
    echo "=========================================="
    echo "Please check the error messages above."
    exit 1
fi

echo ""
echo "=========================================="
echo "Job completed at $(date)"
echo "Total runtime: $SECONDS seconds"
echo "=========================================="