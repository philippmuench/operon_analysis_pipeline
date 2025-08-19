# Operon Analysis Pipeline

A bioinformatics pipeline for analyzing sequence diversity in the fructoselysine/glucoselysine operon across <em>Enterococcus faecalis</em> genomes.

## Overview

This pipeline processes <em>E. faecalis</em> genome assemblies to:
1. Annotate genomes with Prokka
2. Extract reference operon genes from GenBank
3. Search for operon genes across all genomes using BLAST
4. Identify and analyze core genes for comparison
5. Extract operon sequences from assemblies and create MSAs
6. Run diversity analyses and generate plots
7. Run dN/dS and substitution analyses

## Quick Start

### Prerequisites

```bash
# Create conda environment
conda env create -f environment.yml
conda activate efs_diversity
```

### Test Run (50 genomes)

```bash
# 1. Annotate genomes with Prokka
cd 01_prokka_annotation
sbatch run_prokka.sh --test

# 2. Extract operon genes from reference GenBank
cd ../02_reference_operon_extraction
python extract_operon_genes.py

# 3. Search for operon genes via BLAST
cd ../03_blast_search
sbatch run_blast_search.sh
python process_blast_results.py
python create_blast_overview.py
python create_simple_summary.py

# 4. Identify core genes
cd ../04_core_gene_analysis
sbatch run_core_analysis.sh

# 5. Extract operon sequences from assemblies and create MSAs
cd ../05_operon_assembly_extraction
sbatch run_operon_extraction.sh

# 6. Diversity analysis
cd ../06_diversity_analysis
sbatch run_complete_diversity_analysis.sh

# 7. dN/dS and substitution analyses
cd ../07_dnds_analysis
sbatch run_dnds_analysis.sh
```

### Full Run (8,587 genomes)

Remove the `--test` flag from the commands above to process all genomes.

## Pipeline Structure

```
operon_analysis/
├── 01_prokka_annotation/           # Genome annotation with Prokka
├── 02_reference_operon_extraction/ # Extract reference operon genes from GenBank
├── 03_blast_search/                # Search for operon genes via BLAST
├── 04_core_gene_analysis/          # Identify and analyze core genes; create core MSAs
├── 05_operon_assembly_extraction/  # Extract operon sequences from assemblies; create MSAs
├── 06_diversity_analysis/          # Diversity metrics and plots (operon and core)
├── 07_dnds_analysis/               # dN/dS and substitution analyses
├── environment.yml                 # Conda environment specification
└── run_pipeline.sh                 # Master pipeline script
```

## Key Features

- **Modular design**: Each step can be run independently
- **Test mode**: Process 50 genomes for testing
- **SLURM support**: Optimized for HPC clusters
- **Multi-strategy BLAST approaches**: 
  - tblastn (protein→nucleotide) for protein-coding genes
  - blastn (nucleotide→nucleotide) for regulatory elements
  - Multiple extraction strategies for robust comparison

## Target Operon

The fructoselysine/glucoselysine operon contains 7 genes:
- `frpC`: Fructoselysine-6-phosphate deglycase
- `glpC`: Glucoselysine-6-phosphate deglycase
- `ptsD`: PTS system EIID component
- `ptsC`: PTS system EIIC component
- `ptsB`: PTS system EIIB component
- `ptsA`: PTS system EIIA component
- `fruR`: Sigma-54 dependent transcriptional regulator

Plus regulatory elements:
- Promoter region
- Pribnow box (-10 box)

## Results Summary

### Test Run Results (50 genomes)
- **Operon prevalence:** 96% (48/50 genomes)
- **Complete operons:** 48 genomes with all 7 genes
- **Core genes identified:** 1,269 genes (≥95% prevalence)
- **Promoter found in:** 100% of genomes

### Expected Results (8,587 genomes)
- **Operon prevalence:** ~96% of genomes
- **Conservation:** >98% average identity
- **Selection:** Strong purifying (dN/dS < 1)
- **Core genes:** ~1,225 at 95% threshold

## Requirements

- SLURM cluster with `cpu` partition
- Conda environment: `efs_diversity`
- ~500GB storage for intermediate files
- ~200 CPU hours total for full run

## Data Requirements

- Input genomes: `../Efs_assemblies/*.fasta.gz`
- Reference operon: `02_reference_operon_extraction/operon.gb`

## Citation

If you use this pipeline, please cite:
```
[Citation information to be added]
```

## License

[License information to be added]