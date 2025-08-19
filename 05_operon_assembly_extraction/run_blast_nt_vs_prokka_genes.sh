#!/bin/bash
# Run blastn: Prokka CDS (.ffn) of each genome vs reference operon nt (reproduces old qseq-based variants)

set -euo pipefail

PROKKA_DIR="${1:-../01_prokka_annotation/output/prokka_results}"
REF_DB_DIR="${2:-output/reference_db}"
REF_FASTA_SRC="${3:-../02_reference_operon_extraction/output/operon_genes_nt.fasta}"
OUT_DIR="${4:-output/blast_nt_vs_prokka}"

THREADS=${SLURM_CPUS_PER_TASK:-8}

mkdir -p "$OUT_DIR" "$REF_DB_DIR"

REF_FASTA="$REF_DB_DIR/operon_genes_nt.fasta"
cp -f "$REF_FASTA_SRC" "$REF_FASTA"

# Make BLAST DB for reference operon genes (nt)
if [ ! -f "$REF_DB_DIR/operon_genes_nt.nsq" ]; then
  makeblastdb -in "$REF_FASTA" -dbtype nucl -out "$REF_DB_DIR/operon_genes_nt"
fi

mapfile -t GENOME_DIRS < <(find "$PROKKA_DIR" -mindepth 1 -maxdepth 1 -type d | sort)

echo "Running blastn (Prokka .ffn → operon gene DB) for ${#GENOME_DIRS[@]} genomes..."

for GD in "${GENOME_DIRS[@]}"; do
  [ -d "$GD" ] || continue
  FFN=$(ls "$GD"/*.ffn 2>/dev/null | head -n1)
  if [ -z "$FFN" ]; then
    echo "⚠️  Missing .ffn in $GD, skipping"
    continue
  fi
  GENOME=$(basename "$FFN" .ffn)
  OUT_FILE="$OUT_DIR/${GENOME}_blast.txt"
  if [ -s "$OUT_FILE" ]; then
    continue
  fi
  echo "  blastn → $GENOME"
  blastn \
    -query "$FFN" \
    -db "$REF_DB_DIR/operon_genes_nt" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" \
    -perc_identity 80 \
    -num_threads "$THREADS" \
    -out "$OUT_FILE"
done

echo "✅ blastn (Prokka .ffn vs operon DB) completed: $OUT_DIR"


