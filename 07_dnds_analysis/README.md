# Step 7: dN/dS and Substitution Analysis

This directory contains a consolidated pipeline for analyzing selection pressure and substitution patterns in operon and core genes.

## Overview

The consolidated pipeline (`dnds_pipeline.py`) performs comprehensive analysis of:
- **dN/dS ratios** (ratio of non-synonymous to synonymous substitutions)
- **SNP types** (synonymous, non-synonymous)
- **Selection pressure** (purifying, neutral, or positive selection)
- **Sequence variation** (variable sites, parsimony informative sites)
- **Statistical comparisons** between operon and core genes

## Quick Start

```bash
# Run the consolidated pipeline via SLURM
sbatch run_dnds_analysis.sh

# Or run directly with custom parameters
python dnds_pipeline.py \
    --operon-dir ../05_operon_assembly_extraction/output/msa/dna_alignments \
    --core-dir ../04_core_gene_analysis/output/core_gene_alignments \
    --output-dir output \
    --threads 20
```

## Input Data

The pipeline automatically detects and uses MSAs from previous steps:

1. **Operon Gene MSAs** (from Step 05):
   - Primary: `../05_operon_assembly_extraction/output/msa/dna_alignments/`
   - Alternative paths checked automatically

2. **Core Gene MSAs** (from Step 04):
   - Directory: `../04_core_gene_analysis/output/core_gene_alignments/`
   - All available alignments are processed

## Main Components

### Consolidated Pipeline: `dnds_pipeline.py`
Single comprehensive script that:
- Processes all alignments in parallel
- Calculates dN/dS ratios and substitution statistics
- Generates publication-ready visualizations
- Creates detailed reports and summaries
- Performs statistical comparisons

### SLURM Script: `run_dnds_analysis.sh`
Enhanced submission script that:
- Validates input directories
- Checks environment dependencies
- Provides detailed progress reporting
- Shows summary statistics upon completion

## Output Structure

```
output/
├── tables/
│   ├── operon_dnds_analysis.csv      # Detailed operon gene metrics
│   └── core_genes_dnds_analysis.csv  # Detailed core gene metrics
├── plots/
│   ├── dnds_analysis_summary.png     # Comprehensive 6-panel figure
│   └── dnds_analysis_summary.pdf     # PDF version
├── dnds_summary_report.txt           # Comprehensive text report
├── dnds_summary_statistics.txt       # Quick summary statistics
└── dnds_pipeline.log                  # Processing log
```

## Visualizations

The pipeline generates a comprehensive 6-panel figure including:
1. **dN/dS Distribution**: Histogram comparison between gene sets
2. **Substitution Types**: Bar plot of synonymous vs non-synonymous
3. **Selection Categories**: Genes under different selection pressures
4. **Variable Sites**: Correlation with alignment length
5. **Gap Percentages**: Distribution of alignment gaps
6. **Statistical Summary**: Key metrics and test results

## Key Metrics

### dN/dS Ratio Categories
- **< 0.1**: Strong purifying selection
- **0.1-0.5**: Purifying selection
- **0.5-1.0**: Weak purifying selection
- **≈ 1.0**: Neutral evolution
- **> 1.0**: Positive selection

### Additional Metrics
- **Variable sites**: Positions with sequence variation
- **Parsimony informative**: Sites useful for phylogeny
- **Singleton sites**: Unique mutations
- **Gap percentage**: Alignment quality indicator

## Statistical Analysis

The pipeline performs:
- **Mann-Whitney U test**: Non-parametric comparison of dN/dS distributions
- **Effect size**: Rank-biserial correlation
- **Descriptive statistics**: Mean, median, quartiles
- **Selection pressure classification**: Categorization by dN/dS ranges

## Resource Requirements

- **Time**: 2 hours
- **Memory**: 32GB
- **CPUs**: 20
- **Partition**: cpu

## Dependencies

Required Python packages:
- BioPython (≥1.79)
- pandas (≥1.3)
- numpy (≥1.21)
- matplotlib (≥3.4)
- seaborn (≥0.11)
- scipy (≥1.7)

## Troubleshooting

If the pipeline fails:
1. Check input directories exist and contain alignments
2. Review the pipeline log: `output/dnds_pipeline.log`
3. Check SLURM error logs: `dnds_pipeline_*.err`
4. Verify environment activation and package installation
5. Ensure sufficient memory for large datasets