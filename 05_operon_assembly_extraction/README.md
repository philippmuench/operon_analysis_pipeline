# Step 05: Operon Sequence Extraction and Multiple Sequence Alignment

This step extracts operon gene sequences from genome assemblies based on BLAST search results and creates multiple sequence alignments (MSAs) for downstream phylogenetic and evolutionary analyses.

## Overview

The pipeline performs **Strategy D** (assembly-based tblastn) as the primary analysis:
- Extracts sequences directly from raw genome assemblies using BLAST coordinates
- Creates MSAs for both DNA and protein sequences
- Generates conservation plots and metrics
- Optionally includes gene boundary validation using Prokka annotations

## Usage

### Basic Run (Primary Analysis Only)
```bash
sbatch run_operon_extraction.sh
```

### Include Gene Boundary Analysis
```bash
sbatch run_operon_extraction.sh --with-gene-boundary
```

### Run Specific Steps
```bash
# Only create plots and metrics (step 4)
sbatch run_operon_extraction.sh --start-step 4

# Run from MSA creation onwards (steps 2-5)
sbatch run_operon_extraction.sh --start-step 2
```

## Pipeline Steps

1. **Extract sequences from assemblies** - Uses tblastn results to extract DNA sequences directly from genome assemblies
2. **Create MSAs** - Aligns extracted sequences using MAFFT
3. **Extract promoter sequences** - Processes non-coding regulatory regions
4. **Create conservation plots** - Generates Shannon entropy plots and sequence logos
5. **Generate summary** - Creates comprehensive extraction statistics

## Output Structure

```
output/
├── sequences/                    # Extracted DNA sequences (one file per gene)
│   ├── frpC.fasta
│   ├── glpC.fasta
│   ├── ptsA.fasta
│   ├── ptsB.fasta
│   ├── ptsC.fasta
│   ├── ptsD.fasta
│   └── fruR.fasta
├── msa/                          # Multiple sequence alignments
│   ├── dna_alignments/           # DNA alignments for dN/dS analysis
│   ├── protein_alignments/       # Protein alignments
│   └── noncoding_alignments/     # Promoter alignments
├── plots/                        # Conservation visualizations
│   ├── *_shannon_entropy_conservation.png
│   └── *_nucleotide_frequency_logo.png
├── noncoding_sequences/          # Extracted regulatory sequences
├── operon_conservation_metrics.csv  # Conservation statistics
├── extraction_pipeline_summary.txt  # Comprehensive summary
└── manuscript_numbers.txt        # Key statistics for manuscript

# Optional gene boundary analysis
output/gene_boundary_analysis/
├── sequences/                    # Prokka-based extracted sequences
├── msa/                         # Alignments using Prokka boundaries
├── plots/                       # Conservation plots for validation
└── gene_boundary_conservation_metrics.csv
```

## Key Files

### Input Requirements
- BLAST results from step 03: `../03_blast_search/output/blast_results/`
- Genome assemblies: `../../Efs_assemblies/`
- Prokka annotations: `../01_prokka_annotation/output/prokka_results/`

### Main Outputs
- **DNA alignments**: `output/msa/dna_alignments/*_aligned.fasta` - Used for dN/dS analysis
- **Conservation metrics**: `output/operon_conservation_metrics.csv` - Shannon entropy scores
- **Extraction summary**: `output/extraction_pipeline_summary.txt` - Pipeline statistics

## Analysis Strategy

### Primary Analysis (Strategy D)
- **Method**: tblastn (protein → nucleotide)
- **Source**: Raw genome assemblies
- **Purpose**: Phylogenetic and evolutionary analyses
- **Advantages**: 
  - No annotation bias
  - Captures actual genomic sequences
  - Suitable for dN/dS calculations

### Optional Gene Boundary Analysis (Strategy A)
- **Method**: tblastn with Prokka annotations
- **Source**: Prokka-annotated genomes
- **Purpose**: Validate gene boundaries
- **Use**: Quality control only

## Conservation Metrics

The pipeline calculates:
- **Shannon entropy**: Position-specific conservation (0=variable, 1=conserved)
  - For each position in the alignment, counts nucleotide frequencies (A, C, G, T)
  - Calculates Shannon entropy: H = -Σ(p × log₂(p)) where p is the frequency of each nucleotide
  - Converts to conservation score: Conservation = 1 - (H / 2.0) where 2.0 is the maximum entropy for DNA
  - Higher scores indicate more conserved positions
- **Pairwise identity**: Average sequence similarity using efficient count-based calculation
  - For each position, counts occurrences of each nucleotide
  - Calculates number of equal pairs: Σ(count × (count-1) / 2) for each nucleotide
  - Divides by total number of possible pairs: N × (N-1) / 2 where N is number of sequences
  - Averages across all positions to get mean pairwise identity
  - **Performance note**: Uses O(N×L) algorithm instead of O(N²×L) for fast computation with large datasets
- **Gap percentage**: Alignment quality indicator
- **Sequence coverage**: Number of genomes with each gene

## Dependencies

Required tools (available in `efs_diversity` conda environment):
- Python 3.10+ with BioPython
- MAFFT 7.490+
- matplotlib, seaborn for plotting

## Example Commands

```bash
# Full pipeline with gene boundary validation
sbatch run_operon_extraction.sh --with-gene-boundary

# Only regenerate plots from existing alignments
sbatch run_operon_extraction.sh --start-step 4

# Check job status
squeue -u $USER | grep operon

# View results summary
cat output/extraction_pipeline_summary.txt
```

## Pipeline Architecture

### Consolidated Implementation
All extraction and analysis functionality is consolidated in `operon_pipeline.py`:
```bash
# View available commands
python operon_pipeline.py --help

# Run individual steps manually
python operon_pipeline.py extract-sequences --blast-dir ../03_blast_search/output/blast_results \
    --assemblies-dir ../../Efs_assemblies --output-dir output/sequences

python operon_pipeline.py create-msa --sequences-dir output/sequences \
    --output-dir output/msa --threads 8

python operon_pipeline.py create-plots --msa-dir output/msa/dna_alignments \
    --output-dir output/plots

# Require complete operons only (all 7 genes present)
python operon_pipeline.py extract-sequences --require-complete \
    --blast-dir ../03_blast_search/output/blast_results \
    --assemblies-dir ../../Efs_assemblies --output-dir output/sequences
```

### Generate Manuscript Statistics
```bash
# Generate statistics for Methods section
python manuscript_numbers.py
```

## Notes

- The pipeline uses `/vol/tmp` for MAFFT temporary files (faster I/O)
- Extraction thresholds: ≥90% identity, ≥80% coverage
- Promoter extraction uses relaxed thresholds: ≥80% identity, ≥70% coverage
- Processing time: ~4-6 hours for complete pipeline with 8,500+ genomes
- Main script: `run_operon_extraction.sh` (SLURM submission)
- Python implementation: `operon_pipeline.py` (consolidated pipeline with subcommands)
- Statistics generator: `manuscript_numbers.py` (for publication figures)