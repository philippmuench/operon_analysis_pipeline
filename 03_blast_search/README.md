# Step 3: BLAST Search for Operon Genes

This directory contains scripts for searching operon genes across all annotated genomes.

## Input
- **Query**: `../data/operon_genes.fasta` - Protein sequences from the reference operon
- **Databases**: `../prokka_output/*/` - Protein sequences from Prokka annotations

## Output
**Production mode:**
- `../data/blast_results/*_blast.txt` - BLAST results for each genome
- `../results/operon_presence_summary.csv` - Summary of operon presence across genomes
- `../results/all_blast_hits.csv` - Compiled BLAST hits

**Test mode:**
- `output/blast_results/*_blast.txt` - BLAST results for first 50 genomes
- `output/operon_presence_summary.csv` - Summary for test genomes
- `output/all_blast_hits.csv` - Compiled hits for test genomes

## Usage
```bash
# Run BLAST search on all genomes
sbatch run_blast_search.sh

# Test mode: BLAST search on first 50 genomes only
sbatch run_blast_search.sh --test

# Process results (all genomes)
python process_blast_results.py

# Process results (test mode - first 50 genomes)
python process_blast_results.py --test
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