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

### `identify_core_genes.py`
Identifies genes present in ≥95% of genomes from Prokka GFF files.
Automatically detects and reports the number of genomes being processed.

### `extract_core_sequences.py`
Extracts nucleotide sequences for all core genes from Prokka output.
Uses parallel processing to efficiently extract from all genomes.

### `create_core_gene_msa.py`
Creates multiple sequence alignments for all core genes using MAFFT.
Processes genes in parallel for efficiency.

### `run_core_analysis.sh`
SLURM script to run the full analysis pipeline.

## Usage

### Complete Pipeline (Recommended)
```bash
# Run full analysis pipeline via SLURM
sbatch run_core_analysis.sh
```

### Manual Steps
```bash
# Step 1: Identify core genes (processes all genomes in step 01 output)
python identify_core_genes.py

# Step 2: Extract sequences
python extract_core_sequences.py

# Step 3: Create MSAs
python create_core_gene_msa.py
```

## Key Findings (from test run)
- 2,481 unique genes found across 50 test genomes
- 1,269 genes present in ≥95% of test genomes
- 1,058 genes present in 100% of test genomes
- Operon genes (frpC, glpC, etc.) are NOT core genes (present in 96% of test genomes)
## Threshold Analysis

The pipeline now includes comprehensive threshold analysis to understand how different prevalence thresholds affect core gene definitions:

### New Script: plot_core_gene_thresholds.py

**Purpose**: Analyze the impact of different prevalence thresholds on core gene counts

**Features**:
- Creates accumulation curves showing core gene counts vs. prevalence thresholds
- Highlights commonly used thresholds (90%, 95%, 99%, 100%)
- Generates both full-range and zoomed (85-100%) plots
- Produces summary table with gene counts at key thresholds

**Output Files**:
- `core_gene_threshold_curve.png` - Full threshold range (0-100%)
- `core_gene_threshold_curve_zoomed.png` - Zoomed view (85-100%)
- `core_gene_threshold_summary.csv` - Summary table with key statistics

**Usage**:
```bash
python plot_core_gene_thresholds.py \
    --core-genes-file output/gene_prevalence_stats.csv \
    --prokka-dir ../01_prokka_annotation/output/prokka_results \
    --output-dir output
```

**Integration**: Automatically run as part of `run_core_analysis.sh` pipeline.

## Manuscript Statistics

Generate statistics for the manuscript after running core gene analysis:
```bash
# Generate with SLURM (recommended)
sbatch run_manuscript_stats.sh

# Or run directly (outputs to console)
python manuscript_numbers.py

# Save statistics to file
python manuscript_numbers.py manuscript_stats.txt
```

The statistics include:
- Total unique genes identified across genomes
- Number of core genes (≥95% prevalence)
- Conservation scores and categories
- Mean pairwise identity and gap statistics
- Summary ready for manuscript inclusion
