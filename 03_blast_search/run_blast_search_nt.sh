#!/bin/bash
#SBATCH --job-name=blast_operon_nt
#SBATCH --output=output/blast_nt_%A_%a.out
#SBATCH --error=output/blast_nt_%A_%a.err
#SBATCH --array=1-86%20
#SBATCH --time=6:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=32
#SBATCH --partition=cpu

# Nucleotide-vs-nucleotide search for coding genes:
#   Query:  ../02_reference_operon_extraction/output/operon_genes_nt.fasta
#   Subject: Prokka genome .fna
# Output: output/blast_results_nt/<genome>_genes_blast.txt

set -euo pipefail

# Initialize conda (align with run_blast_search.sh)
eval "$(/home/pmuench/miniconda3/bin/conda shell.bash hook)"
conda activate efs_diversity

PROKKA_DIR="../01_prokka_annotation/output/prokka_results"
OUTPUT_DIR="output/blast_results_nt"
QUERY_GENES_NT="../02_reference_operon_extraction/output/operon_genes_nt.fasta"
mkdir -p "$OUTPUT_DIR"

# Batch slicing like the protein script
BATCH_SIZE=100
START_IDX=$(( ($SLURM_ARRAY_TASK_ID - 1) * $BATCH_SIZE + 1 ))
END_IDX=$(( $SLURM_ARRAY_TASK_ID * $BATCH_SIZE ))

mapfile -t ALL_DIRS < <(find "$PROKKA_DIR" -mindepth 1 -maxdepth 1 -type d | sort)
BATCH_DIRS=("${ALL_DIRS[@]:$((START_IDX-1)):$BATCH_SIZE}")

echo "Processing nt batch $SLURM_ARRAY_TASK_ID (${#BATCH_DIRS[@]} genomes)"

detect_prefix() {
  local dir="$1"
  local fna
  fna=$(ls "$dir"/*.fna 2>/dev/null | head -n1)
  if [ -n "$fna" ]; then
    basename "$fna" .fna
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

for GENOME_DIR in "${BATCH_DIRS[@]}"; do
  [ -d "$GENOME_DIR" ] || continue
  GENOME_PREFIX=$(detect_prefix "$GENOME_DIR") || { echo "Skipping $GENOME_DIR (no prefix)"; continue; }
  FNA_FILE="$GENOME_DIR/${GENOME_PREFIX}.fna"
  [ -f "$FNA_FILE" ] || { echo "Missing FNA for $GENOME_PREFIX"; continue; }

  OUT_GENES="$OUTPUT_DIR/${GENOME_PREFIX}_genes_blast.txt"
  if [ -s "$OUT_GENES" ]; then
    echo "Skipping $GENOME_PREFIX (already done)"
    continue
  fi

  echo "blastn nt→genome: $GENOME_PREFIX"
  blastn \
    -query "$QUERY_GENES_NT" \
    -subject "$FNA_FILE" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
    -evalue 1e-5 \
    -num_threads "$SLURM_CPUS_PER_TASK" \
    -out "$OUT_GENES"
done

echo "NT blast completed for batch $SLURM_ARRAY_TASK_ID"



