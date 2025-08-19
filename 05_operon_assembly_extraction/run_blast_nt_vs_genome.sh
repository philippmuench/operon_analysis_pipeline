#!/bin/bash
# Run blastn: reference operon nucleotide CDS vs each genome .fna (Prokka outputs)

set -euo pipefail

PROKKA_DIR="${1:-../01_prokka_annotation/output/prokka_results}"
QUERY_NT="${2:-../02_reference_operon_extraction/output/operon_genes_nt.fasta}"
OUT_DIR="${3:-output/blast_nt_vs_genome}"

THREADS=${SLURM_CPUS_PER_TASK:-8}

mkdir -p "$OUT_DIR"

if [ ! -f "$QUERY_NT" ]; then
  echo "❌ Query file not found: $QUERY_NT" >&2
  exit 1
fi

mapfile -t GENOME_DIRS < <(find "$PROKKA_DIR" -mindepth 1 -maxdepth 1 -type d | sort)

echo "Running blastn (nt vs genome) for ${#GENOME_DIRS[@]} genomes..."

for GD in "${GENOME_DIRS[@]}"; do
  [ -d "$GD" ] || continue
  # Detect prefix (prefer .fna)
  FNA=$(ls "$GD"/*.fna 2>/dev/null | head -n1)
  if [ -z "$FNA" ]; then
    echo "⚠️  Missing .fna in $GD, skipping"
    continue
  fi
  GENOME=$(basename "$FNA" .fna)
  OUT_FILE="$OUT_DIR/${GENOME}_genes_blast.txt"
  if [ -s "$OUT_FILE" ]; then
    continue
  fi
  echo "  blastn → $GENOME"
  blastn \
    -query "$QUERY_NT" \
    -subject "$FNA" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
    -evalue 1e-5 \
    -num_threads "$THREADS" \
    -out "$OUT_FILE"
done

echo "✅ blastn (nt vs genome) completed: $OUT_DIR"


