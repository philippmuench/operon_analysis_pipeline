#!/bin/bash
#SBATCH --job-name=prokka_rerun
#SBATCH --output=prokka_rerun_%j.out
#SBATCH --error=prokka_rerun_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=cpu

# Rerun Prokka for genomes that failed or have incomplete outputs
# Usage: sbatch rerun_failed_genomes.sh [incomplete_genomes_list.txt]

# Load environment
source /vol/projects/BIFO/utils/loadEnv
conda activate efs_diversity

# Get input file
if [ -n "$1" ]; then
    INCOMPLETE_LIST="$1"
else
    # Find most recent incomplete genomes list
    INCOMPLETE_LIST=$(ls -t incomplete_genomes_*.txt 2>/dev/null | head -1)
fi

if [ ! -f "$INCOMPLETE_LIST" ]; then
    echo "Error: No incomplete genomes list found!"
    echo "Usage: sbatch rerun_failed_genomes.sh [incomplete_genomes_list.txt]"
    echo "Or run validate_prokka_outputs.sh first to generate the list"
    exit 1
fi

OUTPUT_DIR="../prokka_output"
FAILED_COUNT=$(wc -l < "$INCOMPLETE_LIST")

echo "Reprocessing failed/incomplete genomes"
echo "====================================="
echo "Input list: $INCOMPLETE_LIST"
echo "Number of genomes to process: $FAILED_COUNT"
echo "Started: $(date)"
echo ""

# Process each genome
SUCCESS_COUNT=0
FAIL_COUNT=0

while read -r GENOME; do
    if [ ! -f "$GENOME" ]; then
        echo "Warning: Genome file not found: $GENOME"
        continue
    fi
    
    BASENAME=$(basename "$GENOME" .fasta.gz)
    OUTDIR="$OUTPUT_DIR/$BASENAME"
    
    echo "Processing $BASENAME..."
    
    # Remove any existing incomplete output
    if [ -d "$OUTDIR" ]; then
        echo "  Removing incomplete output directory..."
        rm -rf "$OUTDIR"
    fi
    
    # Extract genome
    TEMP_FASTA="temp_rerun_${SLURM_JOB_ID}_${BASENAME}.fasta"
    zcat "$GENOME" > "$TEMP_FASTA"
    
    # Run Prokka
    if prokka \
        --outdir "$OUTDIR" \
        --prefix "$BASENAME" \
        --kingdom Bacteria \
        --genus Enterococcus \
        --species faecalis \
        --cpus $SLURM_CPUS_PER_TASK \
        --force \
        --quiet \
        "$TEMP_FASTA"; then
        
        echo "  SUCCESS: $BASENAME"
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    else
        echo "  FAILED: $BASENAME"
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
    
    # Clean up
    rm -f "$TEMP_FASTA"
    echo ""
    
done < "$INCOMPLETE_LIST"

echo "Reprocessing Summary"
echo "==================="
echo "Total genomes processed: $((SUCCESS_COUNT + FAIL_COUNT))"
echo "Successful: $SUCCESS_COUNT"
echo "Failed: $FAIL_COUNT"
echo "Completed: $(date)"

# Run validation check
echo ""
echo "Running validation check..."
cd $(dirname "$0")
./validate_prokka_outputs.sh