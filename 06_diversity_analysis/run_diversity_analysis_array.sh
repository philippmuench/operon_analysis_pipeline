#!/bin/bash
#SBATCH --job-name=diversity_analysis
#SBATCH --output=logs/diversity_%A_%a.out
#SBATCH --error=logs/diversity_%A_%a.err
#SBATCH --array=1-7
#SBATCH --time=2:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2

# SLURM array job to analyze MSA diversity in parallel

# Activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate efs_diversity

# Create output directory
mkdir -p diversity_analysis
mkdir -p logs

# Get list of aligned fasta files
ALIGNED_FILES=(msa_output/*_aligned.fasta)

# Check if any files exist
if [ ${#ALIGNED_FILES[@]} -eq 0 ] || [ ! -f "${ALIGNED_FILES[0]}" ]; then
    echo "Error: No aligned FASTA files found in msa_output/"
    echo "Please run the MSA generation script first (run_msa_and_tree_array.sh)"
    exit 1
fi

# Get the current file based on array task ID (0-indexed)
CURRENT_FILE=${ALIGNED_FILES[$SLURM_ARRAY_TASK_ID-1]}

if [ ! -f "$CURRENT_FILE" ]; then
    echo "Error: Aligned file $CURRENT_FILE not found"
    exit 1
fi

# Extract base filename
BASENAME=$(basename "$CURRENT_FILE" _aligned.fasta)

echo "==============================================="
echo "Processing: $CURRENT_FILE"
echo "Output basename: $BASENAME"
echo "Start time: $(date)"
echo "==============================================="

# Step 1: Run static visualization analysis
echo ""
echo "Step 1: Creating static diversity plots..."
echo "----------------------------------------"

python analyze_msa_diversity.py "$CURRENT_FILE" \
    -o "diversity_analysis/${BASENAME}" \
    -n "$BASENAME"

if [ $? -eq 0 ]; then
    echo "Static analysis completed successfully"
else
    echo "Error: Static analysis failed for $BASENAME"
    exit 1
fi

# Step 2: Run interactive visualization
echo ""
echo "Step 2: Creating interactive HTML visualizations..."
echo "--------------------------------------------------"

python create_interactive_msa_visualization.py "$CURRENT_FILE" \
    -o "diversity_analysis/${BASENAME}_interactive.html" \
    -n "$BASENAME" \
    -w 50

if [ $? -eq 0 ]; then
    echo "Interactive visualization completed successfully"
else
    echo "Error: Interactive visualization failed for $BASENAME"
    exit 1
fi

echo ""
echo "==============================================="
echo "Analysis complete for: $BASENAME"
echo "Output files:"
echo "  - diversity_analysis/${BASENAME}_diversity_analysis.png"
echo "  - diversity_analysis/${BASENAME}_diversity_analysis.pdf"
echo "  - diversity_analysis/${BASENAME}_position_statistics.csv"
echo "  - diversity_analysis/${BASENAME}_summary_report.txt"
echo "  - diversity_analysis/${BASENAME}_interactive.html"
echo "  - diversity_analysis/${BASENAME}_interactive_heatmap.html"
echo "  - diversity_analysis/${BASENAME}_interactive_sliding_window.html"
echo "End time: $(date)"
echo "==============================================="