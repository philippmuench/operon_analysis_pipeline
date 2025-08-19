#!/bin/bash
# Orchestrate multiple extraction/alignment strategies to compare outputs

set -euo pipefail

BASE_DIR=$(cd "$(dirname "$0")" && pwd)
cd "$BASE_DIR"

# Ensure MAFFT temporary directory uses fast local scratch by default
MAFFT_TMP_DEFAULT="/vol/tmp"
mkdir -p "${MAFFT_TMPDIR:-$MAFFT_TMP_DEFAULT}" || true
export MAFFT_TMPDIR="${MAFFT_TMPDIR:-$MAFFT_TMP_DEFAULT}"
export TMPDIR="${TMPDIR:-$MAFFT_TMPDIR}"
echo "Using MAFFT_TMPDIR=$MAFFT_TMPDIR (TMPDIR=$TMPDIR)"

mkdir -p output/strategies
mkdir -p output/mappings/aa_nt_mapping/prokka
mkdir -p output/mappings/aa_nt_mapping/assemblies
mkdir -p output/mappings/nt_nt_mapping/prokka_genome
mkdir -p output/mappings/nt_nt_mapping/prokka_variants

THREADS=${SLURM_CPUS_PER_TASK:-8}

echo "=== Strategy A: tblastn (aa→nt) using Prokka genomes (current) ==="
bash -lc "python extract_operon_sequences.py \
  --prokka_dir ../01_prokka_annotation/output/prokka_results \
  --blast_dir ../03_blast_search/output/blast_results \
  --output_dir output/mappings/aa_nt_mapping/prokka/sequences \
  --min_identity 90 --min_coverage 80 --source prokka"
python create_msa.py --coding-sequences output/mappings/aa_nt_mapping/prokka/sequences --output-dir output/mappings/aa_nt_mapping/prokka/msa --threads "$THREADS"
python create_gene_conservation_plots.py --msa-dir output/mappings/aa_nt_mapping/prokka/msa/dna_alignments --output-dir output/mappings/aa_nt_mapping/prokka/plots --title-suffix "aa_vs_nt; source=prokka"

echo "=== Strategy B: blastn (nt→nt) Prokka genomes → coordinate extraction ==="
bash run_blast_nt_vs_genome.sh ../01_prokka_annotation/output/prokka_results \
  ../02_reference_operon_extraction/output/operon_genes_nt.fasta \
  output/blast_nt_vs_genome
python extract_operon_sequences.py \
  --prokka_dir ../01_prokka_annotation/output/prokka_results \
  --blast_dir output/blast_nt_vs_genome \
  --output_dir output/mappings/nt_nt_mapping/prokka_genome/sequences \
  --min_identity 90 --min_coverage 80 --source prokka
python create_msa.py --coding-sequences output/mappings/nt_nt_mapping/prokka_genome/sequences --output-dir output/mappings/nt_nt_mapping/prokka_genome/msa --threads "$THREADS"
python create_gene_conservation_plots.py --msa-dir output/mappings/nt_nt_mapping/prokka_genome/msa/dna_alignments --output-dir output/mappings/nt_nt_mapping/prokka_genome/plots --title-suffix "nt_vs_nt; source=prokka_genome"

echo "=== Strategy C: blastn (nt→nt) Prokka .ffn vs reference DB (qseq variants, old look) ==="
bash run_blast_nt_vs_prokka_genes.sh ../01_prokka_annotation/output/prokka_results \
  output/reference_db ../02_reference_operon_extraction/output/operon_genes_nt.fasta \
  output/blast_nt_vs_prokka
python create_msa_variants_from_blast.py \
  --blast-dir output/blast_nt_vs_prokka \
  --output-dir output/mappings/nt_nt_mapping/prokka_variants \
  --threads "$THREADS"
python create_gene_conservation_plots.py --msa-dir output/mappings/nt_nt_mapping/prokka_variants/msa_variants --output-dir output/mappings/nt_nt_mapping/prokka_variants/plots --title-suffix "nt_vs_nt; source=prokka_variants"

echo "=== Strategy D: tblastn (aa→nt) with direct assemblies extraction ==="
# Reuse tblastn outputs; extract using --source assemblies for raw assembly files if available
python extract_operon_sequences.py \
  --prokka_dir ../01_prokka_annotation/output/prokka_results \
  --blast_dir ../03_blast_search/output/blast_results \
  --output_dir output/mappings/aa_nt_mapping/assemblies/sequences \
  --min_identity 90 --min_coverage 80 --source assemblies \
  --assemblies_dir ../Efs_assemblies
python create_msa.py --coding-sequences output/mappings/aa_nt_mapping/assemblies/sequences --output-dir output/mappings/aa_nt_mapping/assemblies/msa --threads "$THREADS"
python create_gene_conservation_plots.py --msa-dir output/mappings/aa_nt_mapping/assemblies/msa/dna_alignments --output-dir output/mappings/aa_nt_mapping/assemblies/plots --title-suffix "aa_vs_nt; source=assemblies"

echo "All strategies completed. Outputs organized under output/mappings/{aa_nt_mapping,nt_nt_mapping}/..."


