# Step 3: BLAST Search for Operon Genes

This directory contains scripts for searching operon genes across all annotated genomes.

## Input
- **Query proteins**: `../02_reference_operon_extraction/output/operon_genes_protein.fasta` - Protein sequences from the reference operon
- **Query non-coding**: `../02_reference_operon_extraction/output/operon_noncoding_nt.fasta` - Non-coding sequences from the reference operon
- **Databases**: `../01_prokka_annotation/output/prokka_results/*/` - Genome sequences from Prokka annotations

## Output
- `output/blast_results/*_genes_blast.txt` - BLAST results for protein sequences (tblastn)
- `output/blast_results/*_noncoding_blast.txt` - BLAST results for non-coding sequences (blastn)
- Summary/compiled outputs are written to `../04_core_gene_analysis/output/`:
  - `operon_presence_summary.csv` - Summary of operon presence across genomes
  - `operon_presence_absence_matrix.csv` - Presence/absence matrix
  - `component_prevalence.csv`, `identity_statistics.csv` - Overview stats
  - `all_blast_hits_complete.csv` - All parsed BLAST hits

## Usage
```bash
# Run BLAST search on all genomes (writes per-genome results to output/blast_results)
sbatch run_blast_search.sh

# Process results (writes summaries into ../04_core_gene_analysis/output)
python process_blast_results.py
python create_simple_summary.py
python create_blast_overview.py
```

## Parameters
- E-value threshold: 1e-5
- Identity threshold: 90%
- Coverage threshold: 80%
- Max target sequences: 5

## Output Format
The summary file contains:
- Genome ID
- Presence/absence of each operon gene
- Percent identity for each gene found
- Overall operon completeness (0-1)
- Number of genes found

## Notes
- Prokka directory naming: `run_blast_search.sh` auto-detects each genome's prefix from `.fna`/`.faa`/`.gff` files and does not require a specific suffix (e.g., no `.result` needed).
- Post-processing does not use a test mode; it processes everything present in `output/blast_results` and writes results into `../04_core_gene_analysis/output/` to keep outputs centralized.