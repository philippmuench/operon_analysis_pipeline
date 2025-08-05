#!/bin/bash
#SBATCH --job-name=prokka_all
#SBATCH --output=output/prokka_%A_%a.out
#SBATCH --error=output/prokka_%A_%a.err
#SBATCH --array=1-86%20
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --partition=cpu

# Parse command line arguments
TEST_MODE=false
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --test) TEST_MODE=true ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Initialize conda
eval "$(/home/pmuench/miniconda3/bin/conda shell.bash hook)"
conda activate efs_diversity

# Get list of genome files
GENOME_DIR="../Efs_assemblies"
# In test mode, save to local output folder
if [ "$TEST_MODE" = true ]; then
    OUTPUT_DIR="output/prokka_results"
else
    OUTPUT_DIR="../prokka_output"
fi
mkdir -p $OUTPUT_DIR

# Create list of genomes if it doesn't exist
if [ ! -f genome_list.txt ]; then
    ls $GENOME_DIR/*.fasta.gz > genome_list.txt
fi

# If test mode, create a smaller genome list
if [ "$TEST_MODE" = true ]; then
    echo "TEST MODE: Processing only first 50 genomes"
    head -50 genome_list.txt > genome_list_test.txt
    GENOME_LIST_FILE="genome_list_test.txt"
    # For test mode, only need 1 array job
    if [ $SLURM_ARRAY_TASK_ID -gt 1 ]; then
        echo "Test mode only needs array task 1. Exiting array task $SLURM_ARRAY_TASK_ID"
        exit 0
    fi
else
    GENOME_LIST_FILE="genome_list.txt"
fi

# Calculate which genomes to process in this task
BATCH_SIZE=100
START_IDX=$(( ($SLURM_ARRAY_TASK_ID - 1) * $BATCH_SIZE + 1 ))
END_IDX=$(( $SLURM_ARRAY_TASK_ID * $BATCH_SIZE ))

# Get genomes for this batch
GENOMES=$(sed -n "${START_IDX},${END_IDX}p" $GENOME_LIST_FILE)

echo "Processing batch $SLURM_ARRAY_TASK_ID (genomes $START_IDX to $END_IDX)"

# Process each genome
for GENOME in $GENOMES; do
    BASENAME=$(basename "$GENOME" .fasta.gz)
    OUTDIR="$OUTPUT_DIR/$BASENAME"
    
    # Skip if already processed
    if [ -f "$OUTDIR/$BASENAME.gff" ]; then
        echo "Skipping $BASENAME - already processed"
        continue
    fi
    
    echo "Processing $BASENAME..."
    
    # Extract genome
    zcat "$GENOME" > temp_${SLURM_ARRAY_TASK_ID}_${BASENAME}.fasta
    
    # Run Prokka
    prokka \
        --outdir "$OUTDIR" \
        --prefix "$BASENAME" \
        --kingdom Bacteria \
        --genus Enterococcus \
        --species faecalis \
        --cpus $SLURM_CPUS_PER_TASK \
        --centre X \
        --compliant \
        --force \
        --quiet \
        temp_${SLURM_ARRAY_TASK_ID}_${BASENAME}.fasta
    
    # Clean up
    rm temp_${SLURM_ARRAY_TASK_ID}_${BASENAME}.fasta
done

echo "Batch $SLURM_ARRAY_TASK_ID completed"