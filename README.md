# Operon Conservation Analysis Pipeline

This repository contains the complete analysis pipeline for studying the conservation and diversity of the fructoselysine/glucoselysine utilization operon across 8,587 *Enterococcus faecalis* genomes.

## Pipeline Overview

The analysis consists of individual steps that should be executed sequentially. Each step is self-contained in its own directory with documentation and scripts.

## Analysis Steps

Individual analysis steps are organized in numbered directories from `01_prokka_annotation` to `10_phylophlan`. Each directory contains:
- A README.md file describing the specific analysis
- Pipeline scripts (Python and/or shell scripts)
- Input/output specifications
- Results in the `output/` subdirectory

Please navigate to each directory and execute the steps individually according to their documentation.

## Directory Structure

```
01_prokka_annotation/     - Genome annotation with Prokka
02_reference_operon_extraction/ - Extract reference operon sequences
03_blast_search/          - BLAST-based homology searches
04_core_gene_analysis/    - Core genome analysis for comparison
05_operon_assembly_extraction/ - Extract operon sequences from assemblies
06_diversity_analysis/    - Comparative conservation analysis
07_dnds_analysis/         - Selection pressure (dN/dS) analysis
07b_codon_variation_analysis/ - Codon-level variation and pN/pS analysis
08_start_site_analysis/   - Start codon usage patterns
09_operon_order_analysis/ - Gene order and synteny analysis
10_phylophlan/            - Phylogenetic analysis
```

## Requirements

See `environment.yml` for the complete conda environment specification. Key dependencies include:
- Prokka 1.14.6
- BLAST+ 2.15.0
- MAFFT 7.520
- Python 3.10 with BioPython, pandas, numpy, matplotlib

## Usage

Each step should be executed individually in order:

```bash
cd 01_prokka_annotation
bash run_prokka_pipeline.sh

cd ../02_reference_operon_extraction
python extract_operon_genes.py

# Continue with subsequent steps...
```

## Data

Input genomes and metadata are not included in this repository due to size constraints.

## Results

Key results and visualizations are available in each step's `output/` directory. See `manuscript.md` for integrated results and interpretation.