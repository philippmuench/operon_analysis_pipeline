#!/bin/bash
#SBATCH --job-name=blast_variants
#SBATCH --output=output/blast_variants_%A_%a.out
#SBATCH --error=output/blast_variants_%A_%a.err
#SBATCH --array=1-86%20
#SBATCH --time=3:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpu

# Prokka CDS (.ffn) vs operon nt DB to generate qseq-based variant fragments (old-look)

set -euo pipefail

eval "$(/home/pmuench/miniconda3/bin/conda shell.bash hook)"
conda activate efs_diversity

PROKKA_DIR="../01_prokka_annotation/output/prokka_results"
QUERY_OPERON_NT_SRC="../02_reference_operon_extraction/output/operon_genes_nt.fasta"
REF_DB_DIR="output/reference_db"
OUTPUT_DIR="output/blast_results_prokka_variants"

mkdir -p "$OUTPUT_DIR" "$REF_DB_DIR"

# Prepare reference DB once
cp -f "$QUERY_OPERON_NT_SRC" "$REF_DB_DIR/operon_genes_nt.fasta"
if [ ! -f "$REF_DB_DIR/operon_genes_nt.nsq" ]; then
  makeblastdb -in "$REF_DB_DIR/operon_genes_nt.fasta" -dbtype nucl -out "$REF_DB_DIR/operon_genes_nt"
fi

BATCH_SIZE=100
START_IDX=$(( ($SLURM_ARRAY_TASK_ID - 1) * $BATCH_SIZE + 1 ))
END_IDX=$(( $SLURM_ARRAY_TASK_ID * $BATCH_SIZE ))

mapfile -t ALL_DIRS < <(find "$PROKKA_DIR" -mindepth 1 -maxdepth 1 -type d | sort)
BATCH_DIRS=("${ALL_DIRS[@]:$((START_IDX-1)):$BATCH_SIZE}")

echo "Processing prokka-variants batch $SLURM_ARRAY_TASK_ID (${#BATCH_DIRS[@]} genomes)"

for GD in "${BATCH_DIRS[@]}"; do
  [ -d "$GD" ] || continue
  FFN=$(ls "$GD"/*.ffn 2>/dev/null | head -n1)
  if [ -z "$FFN" ]; then
    echo "Missing .ffn in $GD, skipping"
    continue
  fi
  GENOME=$(basename "$FFN" .ffn)
  OUT_FILE="$OUTPUT_DIR/${GENOME}_blast.txt"
  if [ -s "$OUT_FILE" ]; then
    echo "Skipping $GENOME (already done)"
    continue
  fi
  echo "blastn: $GENOME (.ffn → operon nt DB)"
  blastn \
    -query "$FFN" \
    -db "$REF_DB_DIR/operon_genes_nt" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" \
    -perc_identity 80 \
    -num_threads "$SLURM_CPUS_PER_TASK" \
    -out "$OUT_FILE"
done

echo "Prokka-variants blast completed for batch $SLURM_ARRAY_TASK_ID"



