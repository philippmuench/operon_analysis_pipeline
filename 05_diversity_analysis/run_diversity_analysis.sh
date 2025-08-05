#!/bin/bash
#SBATCH --job-name=diversity_analysis
#SBATCH --output=logs/diversity_analysis_%j.out
#SBATCH --error=logs/diversity_analysis_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2

# SLURM job to analyze all MSA files sequentially

# Activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate efs_diversity

# Create output directories
mkdir -p diversity_analysis
mkdir -p logs

# Check if MSA files exist
if [ ! -d "msa_output" ] || [ -z "$(ls -A msa_output/*_aligned.fasta 2>/dev/null)" ]; then
    echo "Error: No aligned FASTA files found in msa_output/"
    echo "Please run the MSA generation script first."
    exit 1
fi

echo "Starting diversity analysis for all MSA files..."
echo "=============================================="
echo "Start time: $(date)"
echo ""

# Count total files
total_files=$(ls msa_output/*_aligned.fasta 2>/dev/null | wc -l)
current=0

# Process each aligned file
for aligned_file in msa_output/*_aligned.fasta; do
    if [ -f "$aligned_file" ]; then
        current=$((current + 1))
        basename=$(basename "$aligned_file" _aligned.fasta)
        echo ""
        echo "[$current/$total_files] Processing: $basename"
        echo "----------------------------------------"
        
        # Run static visualization analysis
        echo "Creating static plots..."
        python analyze_msa_diversity.py "$aligned_file" \
            -o "diversity_analysis/${basename}" \
            -n "$basename"
        
        if [ $? -ne 0 ]; then
            echo "Warning: Static analysis failed for $basename"
            continue
        fi
        
        # Run interactive visualization
        echo "Creating interactive HTML visualization..."
        python create_interactive_msa_visualization.py "$aligned_file" \
            -o "diversity_analysis/${basename}_interactive.html" \
            -n "$basename" \
            -w 50
        
        if [ $? -ne 0 ]; then
            echo "Warning: Interactive visualization failed for $basename"
            continue
        fi
        
        echo "Completed: $basename"
    fi
done

echo ""
echo "=============================================="
echo "All analyses complete!"
echo "End time: $(date)"
echo ""
echo "Output files are in the 'diversity_analysis/' directory:"
echo "- *_diversity_analysis.png/pdf - Static visualization plots"
echo "- *_position_statistics.csv - Detailed per-position data"
echo "- *_summary_report.txt - Summary statistics"
echo "- *_interactive.html - Interactive visualization"
echo "- *_heatmap.html - Diversity heatmap"
echo "- *_sliding_window.html - Sliding window analysis"
echo ""
echo "To view HTML files on compute node:"
echo "srun --pty python -m http.server 8000"
echo "Then navigate to http://localhost:8000/diversity_analysis/"