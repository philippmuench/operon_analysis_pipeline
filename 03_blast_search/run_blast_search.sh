#!/bin/bash
#SBATCH --job-name=blast_operon
#SBATCH --output=output/blast_%A_%a.out
#SBATCH --error=output/blast_%A_%a.err
#SBATCH --array=1-86%20
#SBATCH --time=6:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=60
#SBATCH --partition=cpu

# Integrated BLAST runner for operon analyses.
# Modes (comma-separated via --modes):
#  - coding_protein   (tblastn: protein vs genome .fna) → output/blast_results
#  - coding_nt        (blastn: operon nt vs genome .fna) → output/blast_results_nt
#  - prokka_variants  (blastn: Prokka .ffn vs operon nt DB; outputs qseq) → output/blast_results_prokka_variants
#  - noncoding        (blastn: promoter nt vs genome .fna) → output/blast_results

set -euo pipefail

# Default: run all modes
MODES="coding_protein,coding_nt,prokka_variants,noncoding"

usage() {
  echo "Usage: sbatch run_blast_search.sh [--modes coding_protein,coding_nt,prokka_variants,noncoding]"
}

while [[ $# -gt 0 ]]; do
  case $1 in
    --modes)
      MODES="$2"; shift 2;;
    -h|--help)
      usage; exit 0;;
    *)
      echo "Unknown arg: $1"; usage; exit 1;;
  esac
done

has_mode() {
  case ",$MODES," in
    *,"$1",*) return 0;;
    *) return 1;;
  esac
}

# Initialize conda
eval "$(/home/pmuench/miniconda3/bin/conda shell.bash hook)"
conda activate efs_diversity

# Setup paths
PROKKA_DIR="../01_prokka_annotation/output/prokka_results"
OUT_PROT_DIR="output/blast_results"
OUT_NT_DIR="output/blast_results_nt"
OUT_VARIANTS_DIR="output/blast_results_prokka_variants"
REF_DB_DIR="output/reference_db"

QUERY_GENES_PROT="../02_reference_operon_extraction/output/operon_genes_protein.fasta"
QUERY_GENES_NT="../02_reference_operon_extraction/output/operon_genes_nt.fasta"
QUERY_NONCODING="../02_reference_operon_extraction/output/operon_noncoding_nt.fasta"

mkdir -p "$OUT_PROT_DIR" "$OUT_NT_DIR" "$OUT_VARIANTS_DIR" "$REF_DB_DIR"

# Prepare operon nt DB if variants mode is used
if has_mode prokka_variants; then
  cp -f "$QUERY_GENES_NT" "$REF_DB_DIR/operon_genes_nt.fasta"
  if [ ! -f "$REF_DB_DIR/operon_genes_nt.nsq" ]; then
    makeblastdb -in "$REF_DB_DIR/operon_genes_nt.fasta" -dbtype nucl -out "$REF_DB_DIR/operon_genes_nt"
  fi
fi

# Get list of genomes for this batch
BATCH_SIZE=100
START_IDX=$(( ($SLURM_ARRAY_TASK_ID - 1) * $BATCH_SIZE + 1 ))
END_IDX=$(( $SLURM_ARRAY_TASK_ID * $BATCH_SIZE ))

mapfile -t ALL_DIRS < <(find "$PROKKA_DIR" -mindepth 1 -maxdepth 1 -type d | sort)
BATCH_DIRS=("${ALL_DIRS[@]:$((START_IDX-1)):$BATCH_SIZE}")

echo "Processing batch $SLURM_ARRAY_TASK_ID (${#BATCH_DIRS[@]} genomes)"
echo "Modes: $MODES"

detect_prefix() {
  local dir="$1"
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

for GENOME_DIR in "${BATCH_DIRS[@]}"; do
  [ -d "$GENOME_DIR" ] || continue
  GENOME_PREFIX=$(detect_prefix "$GENOME_DIR") || { echo "Skipping $GENOME_DIR (no prefix)"; continue; }
  GENOME="$GENOME_PREFIX"

  FNA_FILE="$GENOME_DIR/${GENOME_PREFIX}.fna"
  FAA_FILE="$GENOME_DIR/${GENOME_PREFIX}.faa"
  FFN_FILE=$(ls "$GENOME_DIR"/*.ffn 2>/dev/null | head -n1 || true)

  # coding_protein: tblastn protein vs genome .fna
  if has_mode coding_protein; then
    OUT_GENES="$OUT_PROT_DIR/${GENOME}_genes_blast.txt"
    if [ ! -s "$OUT_GENES" ]; then
      if [ -f "$FNA_FILE" ] && [ -f "$FAA_FILE" ]; then
        echo "tblastn coding_protein → $GENOME"
        tblastn \
          -query "$QUERY_GENES_PROT" \
          -subject "$FNA_FILE" \
          -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
          -evalue 1e-5 \
          -num_threads "$SLURM_CPUS_PER_TASK" \
          -out "$OUT_GENES"
      else
        echo "Warning: Missing .fna or .faa for $GENOME; skipping coding_protein"
      fi
    fi
  fi

  # coding_nt: blastn operon nt vs genome .fna
  if has_mode coding_nt; then
    OUT_GENES_NT="$OUT_NT_DIR/${GENOME}_genes_blast.txt"
    if [ ! -s "$OUT_GENES_NT" ]; then
      if [ -f "$FNA_FILE" ]; then
        echo "blastn coding_nt → $GENOME"
        blastn \
          -query "$QUERY_GENES_NT" \
          -subject "$FNA_FILE" \
          -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
          -evalue 1e-5 \
          -num_threads "$SLURM_CPUS_PER_TASK" \
          -out "$OUT_GENES_NT"
      else
        echo "Warning: Missing .fna for $GENOME; skipping coding_nt"
      fi
    fi
  fi

  # prokka_variants: blastn Prokka .ffn vs operon nt DB (qseq)
  if has_mode prokka_variants; then
    OUT_VARIANTS="$OUT_VARIANTS_DIR/${GENOME}_blast.txt"
    if [ ! -s "$OUT_VARIANTS" ]; then
      if [ -n "$FFN_FILE" ] && [ -f "$FFN_FILE" ]; then
        echo "blastn prokka_variants → $GENOME"
        blastn \
          -query "$FFN_FILE" \
          -db "$REF_DB_DIR/operon_genes_nt" \
          -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" \
          -perc_identity 80 \
          -num_threads "$SLURM_CPUS_PER_TASK" \
          -out "$OUT_VARIANTS"
      else
        echo "Warning: Missing .ffn for $GENOME; skipping prokka_variants"
      fi
    fi
  fi

  # noncoding: blastn promoter nt vs genome .fna
  if has_mode noncoding; then
    OUT_NONCODING="$OUT_PROT_DIR/${GENOME}_noncoding_blast.txt"
    if [ ! -s "$OUT_NONCODING" ]; then
      if [ -f "$FNA_FILE" ]; then
        echo "blastn noncoding → $GENOME"
        blastn \
          -query "$QUERY_NONCODING" \
          -subject "$FNA_FILE" \
          -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
          -evalue 1e-5 \
          -num_threads "$SLURM_CPUS_PER_TASK" \
          -out "$OUT_NONCODING"
      else
        echo "Warning: Missing .fna for $GENOME; skipping noncoding"
      fi
    fi
  fi
done

echo "Batch $SLURM_ARRAY_TASK_ID completed"