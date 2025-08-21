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
- `core_gene_sequences/`: FASTA files for each core gene
- `core_gene_alignments/`: Multiple sequence alignments for each core gene
- `extraction_summary.txt`: Summary of sequence extraction
- `msa_summary.txt`: Summary of MSA creation

## Scripts

### `core_gene_pipeline.py`
Consolidated pipeline script that performs all analysis steps:
1. Identifies genes present in ≥95% of genomes from Prokka GFF files
2. Extracts nucleotide sequences for all core genes
3. Creates multiple sequence alignments using MAFFT
4. Calculates diversity metrics
5. Generates threshold analysis plots

### `run_pipeline.sh`
SLURM submission script to run the pipeline with configurable options.

### `manuscript_numbers.py`
Generates statistics for the manuscript methods section (run after pipeline completes).

## Usage

### Complete Pipeline (Recommended)
```bash
# Run full analysis pipeline via SLURM
sbatch run_pipeline.sh

# Or with custom options
sbatch run_pipeline.sh --threshold 0.99  # Use 99% threshold
sbatch run_pipeline.sh --start-step 3    # Start from MSA creation
```

### Direct Python Execution
```bash
# Run the complete pipeline
python core_gene_pipeline.py

# Or with custom options
python core_gene_pipeline.py --threshold 0.99 --threads 40
python core_gene_pipeline.py --start-step 3  # Start from MSA creation
python core_gene_pipeline.py --help  # See all options
```

## Threshold Analysis

The consolidated pipeline includes comprehensive threshold analysis (Step 5) to understand how different prevalence thresholds affect core gene definitions:

**Features**:
- Creates accumulation curves showing core gene counts vs. prevalence thresholds
- Highlights commonly used thresholds (90%, 95%, 99%, 100%)
- Generates both full-range and zoomed (85-100%) plots
- Produces summary table with gene counts at key thresholds

**Output Files**:
- `core_gene_threshold_curve.png` - Full threshold range (0-100%)
- `core_gene_threshold_summary.csv` - Summary table with key statistics

This analysis is automatically run as Step 5 of the pipeline.

## Manuscript Statistics

Generate statistics for the manuscript after running core gene analysis:
```bash
# Run directly (outputs to console)
python manuscript_numbers.py

# Save statistics to file
python manuscript_numbers.py > manuscript_stats.txt
```

The statistics include:
- Total unique genes identified across genomes
- Number of core genes (≥95% prevalence)
- Conservation scores and categories
- Mean pairwise identity and gap statistics
- Summary ready for manuscript inclusion
