#!/bin/bash
#SBATCH --job-name=blast_operon
#SBATCH --output=output/blast_%A_%a.out
#SBATCH --error=output/blast_%A_%a.err
#SBATCH --array=1-86%20
#SBATCH --time=6:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=60
#SBATCH --partition=cpu

# No command line arguments needed

# Initialize conda
eval "$(/home/pmuench/miniconda3/bin/conda shell.bash hook)"
conda activate efs_diversity

# Setup paths
PROKKA_DIR="../01_prokka_annotation/output/prokka_results"
OUTPUT_DIR="output/blast_results"
QUERY_GENES="../02_operon_extraction/output/operon_genes_protein.fasta"
QUERY_NONCODING="../02_operon_extraction/output/operon_noncoding_nt.fasta"
mkdir -p $OUTPUT_DIR

# Get list of genomes for this batch
BATCH_SIZE=100
START_IDX=$(( ($SLURM_ARRAY_TASK_ID - 1) * $BATCH_SIZE + 1 ))
END_IDX=$(( $SLURM_ARRAY_TASK_ID * $BATCH_SIZE ))

# Get genome directories
GENOME_DIRS=($(ls -d $PROKKA_DIR/ENT_*.result 2>/dev/null | sort))

BATCH_DIRS=("${GENOME_DIRS[@]:$((START_IDX-1)):$BATCH_SIZE}")

echo "Processing batch $SLURM_ARRAY_TASK_ID (${#BATCH_DIRS[@]} genomes)"

# Process each genome
for GENOME_DIR in "${BATCH_DIRS[@]}"; do
    GENOME_FULL=$(basename "$GENOME_DIR")
    GENOME=${GENOME_FULL%.result}  # Remove .result suffix
    FNA_FILE="$GENOME_DIR/${GENOME_FULL}.fna"
    FAA_FILE="$GENOME_DIR/${GENOME_FULL}.faa"
    
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