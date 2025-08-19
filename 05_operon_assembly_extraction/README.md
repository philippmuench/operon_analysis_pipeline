# Step 5: Operon Assembly-based Extraction and MSA Creation

This directory contains scripts for extracting operon gene and promoter sequences directly from assemblies using BLAST-derived coordinates, and for creating multiple sequence alignments for downstream analysis.

## Overview

This step processes BLAST search results to extract actual sequences for operon genes and create MSAs for diversity analysis. It separates the sequence extraction process from the downstream analysis.

## Input Data

- BLAST results from `../03_blast_search/output/`
- Prokka annotations from `../01_prokka_annotation/output/prokka_results/`
 - Reference sequences from `../02_reference_operon_extraction/output/`

## Output Files

All outputs are written to the `output/` directory:
- `sequences/`: FASTA files for each operon gene
- `msa/`: Multiple sequence alignments
- `noncoding_sequences/`: Non-coding sequences (promoter, etc.)
- `msa/noncoding_alignments/`: MSAs for non-coding regions
- `extraction_summary.txt`: Summary of extraction process

### Multi-strategy mapping outputs (refactored)

To make strategy provenance explicit, multi-strategy results are organized under `output/mappings/`:

```
output/mappings/
├── aa_nt_mapping/
│   ├── prokka/            # tblastn (aa→nt), sequences extracted from Prokka .fna
│   │   ├── sequences/
│   │   ├── msa/
│   │   └── plots/
│   └── assemblies/        # tblastn (aa→nt), sequences extracted from raw assemblies
│       ├── sequences/
│       ├── msa/
│       └── plots/
└── nt_nt_mapping/
    ├── prokka_genome/     # blastn (nt→nt), Prokka .fna as subject
    │   ├── sequences/
    │   ├── msa/
    │   └── plots/
    └── prokka_variants/   # blastn (nt→nt), Prokka .ffn vs reference DB (qseq variant style)
        ├── gene_variants/
        ├── msa_variants/
        └── plots/
```

Plots include provenance and sequence counts in their titles, e.g. `… (n=123) — aa_vs_nt; source=prokka`.

## Scripts

### Core Extraction Scripts

#### `extract_operon_sequences.py`
Extracts coding gene sequences for operon genes from Prokka output based on BLAST results.
- Uses BLAST coordinates to locate genes in genomes
- Extracts actual nucleotide sequences
- Creates FASTA files for each gene

#### `extract_noncoding_sequences.py`
Extracts non-coding sequences (promoter regions) from BLAST results.
- Processes non-coding BLAST hits (e.g., `*_noncoding_blast.txt`)
- Extracts promoter sequences from genome files
- Handles coordinate mapping and strand orientation

#### `extract_promoter_from_blast.py`
Specialized promoter extraction from individual BLAST files.
- Processes promoter-specific BLAST results
- Extracts sequences with quality filtering
- Creates promoter-specific FASTA files

#### `extract_sequences_from_blast.py`
Alternative extraction method that works directly from BLAST results.
- Quick extraction without needing original genome files
- Good for initial analysis and validation

### MSA Creation Scripts

#### `create_msa.py`
Create MSAs for coding and non-coding sequences using MAFFT.
 - Handles both coding and non-coding sequences
 - Parallel processing for efficiency
 - Quality control and validation

#### `create_promoter_msa_from_blast.py`
Specialized script for promoter sequence MSAs.
- Extracts promoter sequences from individual BLAST files
- Creates MSAs specifically for non-coding regions
- Generates conservation plots with Pribnow box highlighting

### Promoter Analysis and Visualization

#### `create_promoter_plot.py`
Creates basic promoter conservation plots.
- Generates position-wise conservation profiles
- Basic visualization for promoter analysis

#### `create_promoter_plot_with_pribnow.py`
Enhanced promoter plotting with Pribnow box annotation.
- Highlights the Pribnow box (-10 consensus sequence)
- Shows conservation patterns around regulatory elements
- Publication-ready plots with detailed annotations

## Usage

### SLURM (HPC) submission
Submit the full step via SLURM (defaults to ALL strategies):

```bash
sbatch 05_operon_assembly_extraction/run_operon_extraction.sh
```

Override resources if needed (examples):

```bash
sbatch -c 32 --mem=48G --time=12:00:00 05_operon_assembly_extraction/run_operon_extraction.sh
```

Monitor job and logs:

```bash
squeue -j <JOBID>
tail -f 05_operon_assembly_extraction/operon_extraction_<JOBID>.out
```

### Strategy selection (optional)
By default, all strategies run: `current, nt_vs_genome, prokka_variants, assemblies`.
To restrict strategies:

```bash
sbatch 05_operon_assembly_extraction/run_operon_extraction.sh --strategies current,prokka_variants
```

### Individual Steps
```bash
# Step 1: Extract coding sequences
python extract_operon_sequences.py --prokka-dir ../01_prokka_annotation/output/prokka_results

# Step 2: Extract non-coding sequences
python extract_noncoding_sequences.py

# Step 3: Create MSAs for coding genes
python create_msa.py --sequences-dir output/sequences
```

### Promoter Analysis
```bash
# Extract and analyze promoter sequences
python extract_noncoding_sequences.py
python create_promoter_msa_from_blast.py

# Create promoter conservation plots
python create_promoter_plot_with_pribnow.py --blast-based
```

## Key Features

- **Multiple extraction methods**: Direct from genomes or from BLAST results
- **Comprehensive sequence types**: Coding genes and non-coding regions
- **Quality control**: Identity and coverage thresholds
- **Parallel processing**: Efficient handling of multiple genes/genomes
- **MSA validation**: Automatic quality checks for alignments

## Integration

This step provides input for:
- `../06_diversity_analysis/`: Conservation and diversity analysis
- `../07_dnds_analysis/`: dN/dS ratio calculations

## Notes

- Extraction depends on successful BLAST searches from step 03
- MSA quality depends on sequence similarity and alignment parameters
- Output sequences are used as input for all downstream diversity analyses