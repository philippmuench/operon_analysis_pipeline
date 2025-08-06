# Step 3: BLAST Search for Operon Genes

This directory contains scripts for searching operon genes across all annotated genomes.

## Input
- **Query proteins**: `../02_operon_extraction/output/operon_genes_protein.fasta` - Protein sequences from the reference operon
- **Query non-coding**: `../02_operon_extraction/output/operon_noncoding_nt.fasta` - Non-coding sequences from the reference operon
- **Databases**: `../01_prokka_annotation/output/prokka_results/*/` - Genome sequences from Prokka annotations

## Output
- `output/blast_results/*_genes_blast.txt` - BLAST results for protein sequences (tblastn)
- `output/blast_results/*_noncoding_blast.txt` - BLAST results for non-coding sequences (blastn)
- `output/operon_presence_summary.csv` - Summary of operon presence across genomes
- `output/all_blast_hits.csv` - Compiled BLAST hits

## Usage
```bash
# Run BLAST search on all genomes
sbatch run_blast_search.sh

# Process results
python process_blast_results.py
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