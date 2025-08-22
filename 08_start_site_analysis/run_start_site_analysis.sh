#!/bin/bash
#SBATCH --job-name=start_site_analysis
#SBATCH --output=start_site_%j.out
#SBATCH --error=start_site_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8
#SBATCH --partition=cpu

# Usage:
#   sbatch run_start_site_analysis.sh          # Run full analysis (default)
#   sbatch run_start_site_analysis.sh full     # Run full analysis explicitly
#   sbatch run_start_site_analysis.sh visualize # Only create plots from existing results

set -euo pipefail

# Parse command line arguments
MODE=${1:-"full"}  # full, visualize

# Initialize conda environment
echo "Activating conda environment..."
eval "$(/home/pmuench/miniconda3/bin/conda shell.bash hook)"
conda activate efs_diversity
echo "Conda environment activated: $CONDA_DEFAULT_ENV"

PROKKA_DIR="../01_prokka_annotation/output/prokka_results"
GENE_FASTA="../02_reference_operon_extraction/output/operon_genes_nt.fasta"
OUT_DIR="./output"
METADATA_FILE="../00_annotation/8587_Efs_metadata_ASbarcode.txt"

mkdir -p "$OUT_DIR"

if [ "$MODE" == "visualize" ]; then
    echo "=========================================="
    echo "Running in visualization-only mode"
    echo "Creating plots from existing results..."
    echo "=========================================="
    
    python analyze_start_sites.py \
      --output_dir "$OUT_DIR" \
      --visualize-only
    
    echo "Visualization complete. Plots saved in $OUT_DIR/plots/"
    
elif [ "$MODE" == "full" ]; then
    echo "Running start-site analysis on $PROKKA_DIR"
    echo "Will also perform metadata stratification using $METADATA_FILE"

    python analyze_start_sites.py \
      --prokka_dir "$PROKKA_DIR" \
      --gene_reference_fasta "$GENE_FASTA" \
      --output_dir "$OUT_DIR" \
      --blast_dir "../03_blast_search/output/blast_results_prokka_variants" \
      --metadata "$METADATA_FILE" \
      --max_workers "$SLURM_CPUS_PER_TASK"

    echo "Done. Output in $OUT_DIR"
else
    echo "ERROR: Unknown mode '$MODE'. Use 'full' or 'visualize'"
    exit 1
fi

# Generate manuscript statistics if analysis was successful
if [ -f "$OUT_DIR/start_site_summary.tsv" ]; then
    echo ""
    echo "=========================================="
    echo "Generating manuscript statistics..."
    echo "=========================================="
    
    python manuscript_numbers.py --output "$OUT_DIR/manuscript_stats.txt"
    
    # Display the statistics
    if [ -f "$OUT_DIR/manuscript_stats.txt" ]; then
        echo ""
        echo "Manuscript Statistics Generated:"
        echo "----------------------------------------"
        cat "$OUT_DIR/manuscript_stats.txt"
        echo "----------------------------------------"
    fi
    
    # Check if metadata stratification was performed
    if [ -f "$OUT_DIR/metadata_stratification_summary.txt" ]; then
        echo ""
        echo "Metadata Stratification Summary:"
        echo "----------------------------------------"
        cat "$OUT_DIR/metadata_stratification_summary.txt"
        echo "----------------------------------------"
        echo ""
        echo "Additional stratification files generated:"
        echo "  - start_codon_by_source_niche.tsv"
        echo "  - ptsA_by_source_niche.tsv"
        echo "  - start_codon_by_country.tsv"
        
        if [ -d "$OUT_DIR/plots" ]; then
            echo ""
            echo "Visualization plots generated:"
            echo "  - plots/start_codon_by_source_niche.png/pdf"
            echo "  - plots/ptsA_start_codon_by_niche.png/pdf"
            echo "  - plots/start_codon_by_country.png/pdf"
            echo "  - plots/ttg_usage_heatmap.png/pdf"
        fi
    fi
else
    echo "WARNING: Analysis results not found, skipping manuscript statistics"
fi


