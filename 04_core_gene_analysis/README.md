# Step 4: Core Gene Analysis

This directory contains scripts for identifying and analyzing core genes across E. faecalis genomes.

## Overview

This step identifies genes that are present in ≥95% of E. faecalis genomes (core genes) to provide a baseline for comparing the diversity of the operon genes.

## Input Data

Uses Prokka annotations from `../01_prokka_annotation/output/prokka_results/`
- If step 01 was run with `--test`: 50 genomes
- If step 01 was run without `--test`: all 8,587 genomes

## Output Files

All outputs are written to the `output/` directory:
- `core_genes_95pct.txt`: List of genes present in ≥95% of genomes
- `gene_prevalence_stats.csv`: Prevalence statistics for all genes
- `core_gene_sequences/`: FASTA files for each core gene (created by extract_core_sequences.sh)

## Scripts

### `identify_core_genes.py`
Identifies genes present in ≥95% of genomes from Prokka GFF files.
Automatically detects and reports the number of genomes being processed.

### `extract_core_sequences.sh` 
Extracts nucleotide sequences for all core genes from each genome.

### `calculate_diversity.py`
Calculates pairwise identity and diversity metrics for core genes.

### `run_core_analysis.sh`
SLURM script to run the full analysis pipeline.

## Usage
```bash
# Step 1: Identify core genes (processes all genomes in step 01 output)
python identify_core_genes.py

# Step 2: Extract sequences
sbatch extract_core_sequences.sh

# Step 3: Calculate diversity
python calculate_diversity.py
```

## Key Findings (from test run)
- 2,481 unique genes found across 50 test genomes
- 1,269 genes present in ≥95% of test genomes
- 1,058 genes present in 100% of test genomes
- Operon genes (frpC, glpC, etc.) are NOT core genes (present in 96% of test genomes)