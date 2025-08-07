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
QUERY_GENES="../02_reference_operon_extraction/output/operon_genes_protein.fasta"
QUERY_NONCODING="../02_reference_operon_extraction/output/operon_noncoding_nt.fasta"
mkdir -p "$OUTPUT_DIR"

# Get list of genomes for this batch
BATCH_SIZE=100
START_IDX=$(( ($SLURM_ARRAY_TASK_ID - 1) * $BATCH_SIZE + 1 ))
END_IDX=$(( $SLURM_ARRAY_TASK_ID * $BATCH_SIZE ))

### Determine genome directories robustly (no reliance on ".result" suffix)
mapfile -t ALL_DIRS < <(find "$PROKKA_DIR" -mindepth 1 -maxdepth 1 -type d | sort)

# Slice the batch
BATCH_DIRS=("${ALL_DIRS[@]:$((START_IDX-1)):$BATCH_SIZE}")

echo "Processing batch $SLURM_ARRAY_TASK_ID (${#BATCH_DIRS[@]} genomes)"

# Helper to detect prefix inside a Prokka output directory
detect_prefix() {
    local dir="$1"
    # Prefer .fna, then .faa, then .gff
    local fna
    fna=$(ls "$dir"/*.fna 2>/dev/null | head -n1)
    if [ -n "$fna" ]; then
        basename "$fna" .fna
        return 0
    fi
    local faa
    faa=$(ls "$dir"/*.faa 2>/dev/null | head -n1)
    if [ -n "$faa" ]; then
        basename "$faa" .faa
        return 0
    fi
    local gff
    gff=$(ls "$dir"/*.gff 2>/dev/null | head -n1)
    if [ -n "$gff" ]; then
        basename "$gff" .gff
        return 0
    fi
    return 1
}

# Process each genome
for GENOME_DIR in "${BATCH_DIRS[@]}"; do
    [ -d "$GENOME_DIR" ] || continue
    GENOME_PREFIX=$(detect_prefix "$GENOME_DIR")
    if [ -z "$GENOME_PREFIX" ]; then
        echo "Warning: Could not detect Prokka prefix in $GENOME_DIR (missing .fna/.faa/.gff). Skipping."
        continue
    fi

    FNA_FILE="$GENOME_DIR/${GENOME_PREFIX}.fna"
    FAA_FILE="$GENOME_DIR/${GENOME_PREFIX}.faa"

    GENOME="$GENOME_PREFIX"

    if [ ! -f "$FNA_FILE" ] || [ ! -f "$FAA_FILE" ]; then
        echo "Warning: Missing .fna or .faa for $GENOME in $GENOME_DIR"
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
    if [ ! -f "$OUTPUT_GENES" ]; then
        tblastn \
            -query "$QUERY_GENES" \
            -subject "$FNA_FILE" \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
            -evalue 1e-5 \
            -num_threads "$SLURM_CPUS_PER_TASK" \
            -out "$OUTPUT_GENES"
    fi

    # Run BLASTN for non-coding sequences (nucleotide vs nucleotide)
    if [ ! -f "$OUTPUT_NONCODING" ]; then
        blastn \
            -query "$QUERY_NONCODING" \
            -subject "$FNA_FILE" \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
            -evalue 1e-5 \
            -num_threads "$SLURM_CPUS_PER_TASK" \
            -out "$OUTPUT_NONCODING"
    fi
done

echo "Batch $SLURM_ARRAY_TASK_ID completed"