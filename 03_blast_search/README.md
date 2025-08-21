# Step 3: BLAST Search for Operon Genes

This directory contains a unified pipeline for searching operon genes across all annotated genomes using multiple BLAST strategies.

## Pipeline Overview

The unified `blast_pipeline.py` combines BLAST search, result processing, and statistics generation into a single, efficient workflow. It's executed through `run_pipeline.sh` which handles SLURM job submission and mode selection.

## Input
- **Query proteins**: `../02_reference_operon_extraction/output/operon_genes_protein.fasta` (for tblastn)
- **Query coding (nt)**: `../02_reference_operon_extraction/output/operon_genes_nt.fasta` (for blastn coding)
- **Query non-coding**: `../02_reference_operon_extraction/output/operon_noncoding_nt.fasta` (for promoter blastn)
- **Databases**: 
  - `../../Efs_assemblies/*.fasta.gz` - Genome assemblies (main BLAST searches)
  - `../01_prokka_annotation/output/prokka_results/*/` - Prokka annotations (only for variant analysis)

## Output
- `output/blast_results/*_genes_blast.txt` - BLAST results for coding genes (tblastn against assemblies)
- `output/blast_results_nt/*_genes_blast.txt` - BLAST results for coding genes (blastn against assemblies)
- `output/blast_results_prokka_variants/*_blast.txt` - BLAST results for Prokka gene variants (blastn of Prokka genes vs reference)
- `output/blast_results/*_noncoding_blast.txt` - BLAST results for non-coding sequences (blastn against assemblies)
- Summary/compiled outputs are written to `../04_core_gene_analysis/output/`:
  - `operon_presence_summary.csv` - Summary of operon presence across genomes
  - `operon_presence_absence_matrix.csv` - Presence/absence matrix
  - `component_prevalence.csv`, `identity_statistics.csv` - Overview stats
  - `all_blast_hits_complete.csv` - All parsed BLAST hits

## Usage

### Running the Complete Pipeline

```bash
# Step 1: Run BLAST searches across all genomes (array job, processes 100 genomes per batch)
sbatch --array=1-86 run_pipeline.sh

# Step 2: After all array jobs complete, process results
sbatch --array=1 run_pipeline.sh process

# Step 3: Generate overview statistics
sbatch --array=1 run_pipeline.sh overview

# Step 4: Generate manuscript statistics
sbatch --array=1 run_pipeline.sh stats

# Alternative: Run steps 2-4 together after BLAST search completes
sbatch --array=1 run_pipeline.sh all
```

### Pipeline Modes

- **search**: Run BLAST searches (requires array job)
- **process**: Process BLAST results and create summaries
- **overview**: Create comprehensive BLAST overview
- **stats**: Generate manuscript statistics
- **all**: Run complete analysis (process + overview + stats)
- **test**: Test mode with subset of genomes

### Direct Python Execution (Advanced)

```bash
# Run BLAST search for specific batch
python blast_pipeline.py --mode search --batch-id 1 --batch-size 100 --threads 60

# Process all results
python blast_pipeline.py --mode process --threads 60

# Generate statistics
python blast_pipeline.py --mode stats
```

## BLAST Parameters
- E-value threshold: 1e-5
- Identity threshold: 90%
- Coverage threshold: 80%
- Max target sequences: 5
- BLAST modes: coding_protein, coding_nt, prokka_variants, noncoding

## Performance Configuration
- Batch size: 100 genomes per array job
- CPUs per task: 60
- Memory: 128G per job
- Time limit: 6 hours per batch
- Max parallel jobs: 20 (controlled by --array=1-86%20)

## Output Format
The summary files contain:
- Genome ID
- Presence/absence of each operon gene
- Percent identity for each gene found
- Overall operon completeness (0-1)
- Number of genes found
- BLAST hit details (e-value, coverage, alignment positions)

## Manuscript Statistics
The statistics generated include:
- Number of BLAST hits across all search strategies
- Hit quality metrics (identity, coverage distributions)
- Operon completeness across genomes
- Summary of high-quality hits
- Comparative analysis between different BLAST modes

## Notes
- The pipeline automatically detects genome prefixes from `.fna`/`.faa`/`.gff` files
- Array jobs process genomes in batches for efficiency
- All outputs are centralized in `../04_core_gene_analysis/output/` for downstream analysis
- The pipeline includes comprehensive error handling and logging
- Progress is tracked in `output/pipeline_*.out` and `output/pipeline_*.err` files