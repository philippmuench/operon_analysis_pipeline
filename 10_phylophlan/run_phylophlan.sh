#!/bin/bash
#SBATCH --job-name=phylophlan_isolates
#SBATCH --output=logs/phylophlan_%j.out
#SBATCH --error=logs/phylophlan_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=80
#SBATCH --partition=cpu

echo "=========================================="
echo "PhyloPhlAn Analysis for E. faecalis Isolates"
echo "Job ID: $SLURM_JOB_ID"
echo "Date: $(date)"
echo "Cores: $SLURM_CPUS_PER_TASK"
echo "Memory: 128G"
echo "=========================================="

# Initialize conda
echo "Activating conda environment..."
eval "$(/home/pmuench/miniconda3/bin/conda shell.bash hook)"
conda activate phylophlan3
echo "Conda environment activated: $CONDA_DEFAULT_ENV"

# Show working directory
echo "Working directory: $(pwd)"
echo ""

# Create logs directory if it doesn't exist
mkdir -p logs

# Check for input directory
echo "Checking for input files..."
if [ -d "input_isolates" ]; then
    FILE_COUNT=$(find input_isolates -name "*.faa" -o -name "*.fna" | wc -l)
    echo "Found $FILE_COUNT input files in input_isolates/"
else
    echo "ERROR: Input directory 'input_isolates' not found"
    echo "Please create the directory and add genome files first"
    exit 1
fi

# Check for config file
if [ ! -f "isolates_config.cfg" ]; then
    echo "ERROR: Configuration file 'isolates_config.cfg' not found"
    echo "Please create the configuration file first"
    exit 1
fi

echo ""
echo "=========================================="
echo "Running PhyloPhlAn Analysis"
echo "=========================================="

# Run PhyloPhlAn with specified parameters
phylophlan \
    -i input_isolates \
    -o output_isolates \
    -d s__Enterococcus_faecalis \
    --trim greedy \
    --not_variant_threshold 0.99 \
    --remove_fragmentary_entries \
    --fragmentary_threshold 0.67 \
    --min_num_entries 50 \
    -t a \
    -f isolates_config.cfg \
    --diversity low \
    --force_nucleotides \
    --nproc $SLURM_CPUS_PER_TASK \
    --verbose 2>&1 | tee logs/phylophlan__output_isolates.log

# Check if PhyloPhlAn succeeded
if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "PhyloPhlAn completed successfully!"
    echo "=========================================="
    
    # Show output summary
    OUTPUT_DIR="output_isolates"
    if [ -d "$OUTPUT_DIR" ]; then
        echo ""
        echo "Generated files:"
        ls -lh "$OUTPUT_DIR"/ 2>/dev/null | head -20
        
        # Check for key output files
        if [ -f "$OUTPUT_DIR/RAxML_bestTree.output_isolates_refined.tre" ]; then
            echo ""
            echo "Phylogenetic tree generated: $OUTPUT_DIR/RAxML_bestTree.output_isolates_refined.tre"
        fi
        
        if [ -f "$OUTPUT_DIR/output_isolates.tree.nwk" ]; then
            echo ""
            echo "Newick tree file: $OUTPUT_DIR/output_isolates.tree.nwk"
        fi
    fi
else
    echo ""
    echo "ERROR: PhyloPhlAn failed!"
    echo "Check the error messages above and the log file at logs/phylophlan__output_isolates.log"
    exit 1
fi

echo ""
echo "=========================================="
echo "Analysis complete at $(date)"
echo "Total runtime: $SECONDS seconds"
echo "=========================================="