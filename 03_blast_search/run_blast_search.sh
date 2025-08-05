#!/bin/bash
#SBATCH --job-name=blast_operon
#SBATCH --output=output/blast_%A_%a.out
#SBATCH --error=output/blast_%A_%a.err
#SBATCH --array=1-87%20
#SBATCH --time=6:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
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

# Setup paths
# In test mode, look for prokka results in test location
if [ "$TEST_MODE" = true ]; then
    PROKKA_DIR="../01_prokka_annotation/output/prokka_results"
    OUTPUT_DIR="output/blast_results"
    QUERY_GENES="../02_operon_extraction/output/operon_genes_protein.fasta"
    QUERY_NONCODING="../02_operon_extraction/output/operon_noncoding_nt.fasta"
else
    PROKKA_DIR="../prokka_output"
    OUTPUT_DIR="../data/blast_results"
    QUERY_GENES="../data/operon_genes_protein.fasta"
    QUERY_NONCODING="../data/operon_noncoding_nt.fasta"
fi
mkdir -p $OUTPUT_DIR

# Get list of genomes for this batch
BATCH_SIZE=100
START_IDX=$(( ($SLURM_ARRAY_TASK_ID - 1) * $BATCH_SIZE + 1 ))
END_IDX=$(( $SLURM_ARRAY_TASK_ID * $BATCH_SIZE ))

# Get genome directories
GENOME_DIRS=($(ls -d $PROKKA_DIR/ENT_* 2>/dev/null | sort))

# If test mode, limit to first 50 genomes
if [ "$TEST_MODE" = true ]; then
    echo "TEST MODE: Processing only first 50 genomes"
    GENOME_DIRS=("${GENOME_DIRS[@]:0:50}")
    # For test mode, only need array task 1
    if [ $SLURM_ARRAY_TASK_ID -gt 1 ]; then
        echo "Test mode only needs array task 1. Exiting array task $SLURM_ARRAY_TASK_ID"
        exit 0
    fi
fi

BATCH_DIRS=("${GENOME_DIRS[@]:$((START_IDX-1)):$BATCH_SIZE}")

echo "Processing batch $SLURM_ARRAY_TASK_ID (${#BATCH_DIRS[@]} genomes)"

# Process each genome
for GENOME_DIR in "${BATCH_DIRS[@]}"; do
    GENOME=$(basename "$GENOME_DIR")
    FNA_FILE="$GENOME_DIR/${GENOME}.fna"
    FAA_FILE="$GENOME_DIR/${GENOME}.faa"
    
    if [ ! -f "$FNA_FILE" ] || [ ! -f "$FAA_FILE" ]; then
        echo "Warning: Missing .fna or .faa file for $GENOME"
        continue
    fi
    
    OUTPUT_GENES="$OUTPUT_DIR/${GENOME}_genes_blast.txt"
    OUTPUT_NONCODING="$OUTPUT_DIR/${GENOME}_noncoding_blast.txt"
    
    # Skip if already processed
    if [ -f "$OUTPUT_GENES" ] && [ -f "$OUTPUT_NONCODING" ]; then
        echo "Skipping $GENOME - already processed"
        continue
    fi
    
    # Run tblastn for gene sequences (protein query vs translated nucleotide database)
    # This searches protein sequences against translated nucleotide database
    if [ ! -f "$OUTPUT_GENES" ]; then
        tblastn \
            -query $QUERY_GENES \
            -subject $FNA_FILE \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
            -evalue 1e-5 \
            -num_threads $SLURM_CPUS_PER_TASK \
            -out $OUTPUT_GENES
    fi
    
    # Run BLASTN for non-coding sequences (nucleotide vs nucleotide)
    if [ ! -f "$OUTPUT_NONCODING" ]; then
        blastn \
            -query $QUERY_NONCODING \
            -subject $FNA_FILE \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
            -evalue 1e-5 \
            -num_threads $SLURM_CPUS_PER_TASK \
            -out $OUTPUT_NONCODING
    fi
done

echo "Batch $SLURM_ARRAY_TASK_ID completed"