# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.
## Github

when commiting, dont indicate that cluade code worked on this project, so just commit with a normal description

## Project Overview

This bioinformatics pipeline analyzes sequence diversity in the fructoselysine/glucoselysine operon across 8,587 Enterococcus faecalis genomes. The project is organized as a 6-step analysis pipeline that processes genome assemblies from annotation through diversity analysis and visualization.

## Key Objectives

1. **Genome Annotation**: Process all E. faecalis genomes with Prokka
2. **Operon Gene Extraction**: Extract reference operon genes from GenBank file
3. **BLAST Search**: Identify operon genes across all genomes
4. **Core Gene Analysis**: Identify and analyze core genes (≥95% prevalence)
5. **Diversity Analysis**: Calculate conservation metrics and dN/dS ratios
6. **Visualization**: Generate publication-ready figures

## Project Structure

```
operon_analysis/
├── 01_prokka_annotation/    # Genome annotation with Prokka
├── 02_operon_extraction/    # Extract reference genes from GenBank
├── 03_blast_search/         # Search for operon genes via BLAST
├── 04_core_gene_analysis/   # Identify and analyze core genes
├── 05_diversity_analysis/   # Calculate diversity metrics
├── 06_visualization/        # Create figures
├── run_pipeline.sh          # Master pipeline script
└── ../Efs_assemblies/       # 8,587 E. faecalis genomes (.fasta.gz)
```

## Environment Setup

**IMPORTANT: We are on the login node - DO NOT run compute-intensive jobs directly!**
- Use interactive sessions with `srun` for testing
- Submit SLURM jobs with `sbatch` for production runs
- Use `--partition=cpu` and allow sufficient time (e.g., `--time=06:00:00`)

```bash
source /vol/projects/BIFO/utils/loadEnv
conda activate efs_diversity
```

The environment includes:
- Prokka 1.14.6 (genome annotation)
- HMMER 3.4 (HMM-based searches)
- MAFFT 7.520 & MUSCLE 5.1 (MSA tools)
- FastTree 2.1.11 (phylogenetic trees)
- Python 3.10 with BioPython, pandas, matplotlib, seaborn
- BLAST+ (via system)

## Common Commands

**Run complete pipeline (72 hours):**
```bash
./run_pipeline.sh
```

**Run individual steps:**
```bash
# Step 1: Prokka annotation (24 hours, array job)
cd 01_prokka_annotation
sbatch run_prokka.sh
./check_progress.sh  # Monitor progress

# Step 2: Extract operon genes (<1 minute)
cd 02_operon_extraction
python extract_operon_genes.py

# Step 3: BLAST search (4 hours)
cd 03_blast_search
sbatch run_blast_search.sh
python process_blast_results.py

# Step 4: Core gene analysis (12 hours)
cd 04_core_gene_analysis
python identify_core_genes.py
sbatch run_core_analysis.sh

# Step 5: Diversity analysis (6 hours)
cd 05_diversity_analysis
sbatch run_diversity_analysis.sh

# Step 6: Generate figures (30 minutes)
cd 06_visualization
python create_all_figures.py
```

**Interactive testing:**
```bash
srun --cpus-per-task=4 --mem=8G --time=06:00:00 --partition=cpu --pty bash
```

## Target Operon Genes

The fructoselysine/glucoselysine operon contains 7 genes:
- `frpC`: Fructoselysine-6-phosphate deglycase
- `glpC`: Glucoselysine-6-phosphate deglycase  
- `ptsD`: PTS system EIID component
- `ptsC`: PTS system EIIC component
- `ptsB`: PTS system EIIB component
- `ptsA`: PTS system EIIA component
- `fruR`: Sigma-54 dependent transcriptional regulator

## SLURM Configuration

- Array jobs process genomes in batches of 100
- Maximum 20 concurrent array tasks (`--array=1-86%20`)
- Standard resources: 4 CPUs, 8GB RAM per task
- Partition: `cpu`

## Key Files and Outputs

**Input Data:**
- `../Efs_assemblies/*.fasta.gz`: Compressed genome assemblies
- `../data/operon.gb`: Reference operon GenBank file

**Intermediate Files:**
- `../prokka_output/*/`: Prokka annotations (GFF, FAA, FFN)
- `../data/blast_results/*_blast.txt`: BLAST results
- `../data/core_genes_95pct.txt`: List of core genes

**Final Results:**
- `../results/operon_presence_summary.csv`: Operon prevalence data
- `../results/core_gene_diversity.csv`: Core gene diversity metrics
- `../results/substitution_analysis.csv`: dN/dS analysis
- `../results/figures/`: Publication figures

## Important Notes

- Genome files are compressed (`.fasta.gz`) - handle decompression in scripts
- Use array jobs for parallel processing of 8,587 genomes
- Check for existing outputs before reprocessing to save time
- BLAST searches use protein sequences (`.faa` files from Prokka)
- Core gene threshold is 95% prevalence across all genomes