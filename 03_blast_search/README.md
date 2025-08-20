# Step 3: BLAST Search for Operon Genes

This directory contains scripts for searching operon genes across all annotated genomes.

## Input
- **Query proteins**: `../02_reference_operon_extraction/output/operon_genes_protein.fasta` (for tblastn)
- **Query coding (nt)**: `../02_reference_operon_extraction/output/operon_genes_nt.fasta` (for blastn coding)
- **Query non-coding**: `../02_reference_operon_extraction/output/operon_noncoding_nt.fasta` (for promoter blastn)
- **Databases**: `../01_prokka_annotation/output/prokka_results/*/` - Genome sequences from Prokka annotations

## Output
- `output/blast_results/*_genes_blast.txt` - BLAST results for coding genes (tblastn)
- `output/blast_results_nt/*_genes_blast.txt` - BLAST results for coding genes (blastn, nt vs genome)
- `output/blast_results/*_noncoding_blast.txt` - BLAST results for non-coding sequences (blastn)
- Summary/compiled outputs are written to `../04_core_gene_analysis/output/`:
  - `operon_presence_summary.csv` - Summary of operon presence across genomes
  - `operon_presence_absence_matrix.csv` - Presence/absence matrix
  - `component_prevalence.csv`, `identity_statistics.csv` - Overview stats
  - `all_blast_hits_complete.csv` - All parsed BLAST hits

## Usage
```bash
# Default: run ALL modes (coding_protein, coding_nt, prokka_variants, noncoding)
sbatch run_blast_search.sh

# Optionally choose modes explicitly (comma-separated)
sbatch run_blast_search.sh --modes coding_protein,coding_nt

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

## Manuscript Statistics
Generate statistics for the manuscript after running BLAST searches:
```bash
# Generate with SLURM (recommended)
sbatch run_manuscript_stats.sh

# Or run directly (outputs to console)
python manuscript_numbers.py

# Save statistics to file
python manuscript_numbers.py manuscript_stats.txt
```

The statistics include:
- Number of BLAST hits across all search strategies
- Hit quality metrics (identity, coverage)
- Operon completeness across genomes
- Summary of high-quality hits

## Notes
- Prokka directory naming: `run_blast_search.sh` auto-detects each genome's prefix from `.fna`/`.faa`/`.gff` files and does not require a specific suffix (e.g., no `.result` needed).
- Post-processing does not use a test mode; it processes everything present in `output/blast_results` and writes results into `../04_core_gene_analysis/output/` to keep outputs centralized.