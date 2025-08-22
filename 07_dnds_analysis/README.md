# Step 7: dN/dS and Substitution Analysis

This directory contains an optimized pipeline for analyzing selection pressure and substitution patterns in operon and core genes using the NG86 method with proper site normalization.

## Overview

The pipeline (`dnds_pipeline.py`) performs comprehensive analysis of:
- **dN/dS ratios** using NG86 method with site normalization and Jukes-Cantor correction
- **Reference-based analysis** for large alignments (>100 sequences)
- **Pairwise analysis** for smaller alignments
- **Selection pressure classification** (purifying, neutral, or positive selection)
- **Statistical comparisons** between operon and core genes

## Key Features

### Optimized NG86 Implementation
- Proper synonymous/non-synonymous site normalization
- Jukes-Cantor correction for multiple hits
- LRU caching (4096 entries) for efficient computation
- Reference-based approach using actual reference sequences from Step 02

### Site-Based SNP Counting (NEW)
- Simple, intuitive codon-level variation metrics
- Counts variable codon positions as synonymous or non-synonymous
- No complex models - direct interpretation for main text figures
- Complements NG86 for comprehensive analysis

### Performance Optimizations
- Automatic method selection based on alignment size
- Aggregated rate calculations instead of individual omega means
- Optional simple method computation with `--include-simple` flag
- Efficient handling of large alignments with thousands of sequences
- Fixed handling of gap-heavy core gene alignments

## Quick Start

```bash
# Analyze both operon and core genes (default)
sbatch run_dnds_analysis.sh

# Analyze only operon genes
sbatch run_dnds_analysis.sh --mode operon

# Analyze only core genes  
sbatch run_dnds_analysis.sh --mode core

# Include simple method comparison (slower)
sbatch run_dnds_analysis.sh --include-simple

# Custom parameters via direct execution
python dnds_pipeline.py \
    --mode operon \
    --operon-dir ../05_operon_assembly_extraction/output/msa/dna_alignments \
    --output-dir output \
    --threads 20
```

## Command-Line Options

### Analysis Mode
- `--mode {both,operon,core}`: Choose which gene sets to analyze
  - `both` (default): Analyze both operon and core genes
  - `operon`: Analyze only operon genes
  - `core`: Analyze only core genes

### Input/Output
- `--operon-dir`: Directory containing operon gene alignments
- `--core-dir`: Directory containing core gene alignments  
- `--output-dir`: Output directory for results (default: `output`)

### Performance
- `--threads`: Number of threads for parallel processing (default: 20)
- `--include-simple`: Also compute naive count ratio for comparison (slower)

## Input Data

The pipeline automatically detects and uses MSAs from previous steps:

1. **Operon Gene MSAs** (from Step 05):
   - Directory: `../05_operon_assembly_extraction/output/msa/dna_alignments/`
   - Files: `*_aligned.fasta` (e.g., `frpC_aligned.fasta`)

2. **Core Gene MSAs** (from Step 04):
   - Directory: `../04_core_gene_analysis/output/core_gene_alignments/`
   - All available alignments are processed

3. **Reference Sequences** (from Step 02):
   - File: `../02_reference_operon_extraction/output/operon_genes_nt.fasta`
   - Used for reference-based dN/dS calculation

## Output Structure

```
output/
├── tables/
│   ├── operon_dnds_analysis.csv      # Detailed operon gene metrics
│   └── core_genes_dnds_analysis.csv  # Detailed core gene metrics
├── plots/
│   ├── dnds_analysis_summary.png     # Comprehensive 6-panel figure
│   ├── dnds_analysis_summary.pdf     # PDF version
│   └── method_comparison.png         # Comparison of methods (if --include-simple)
├── dnds_summary_report.txt           # Comprehensive text report
├── dnds_summary_statistics.txt       # Quick summary statistics
└── dnds_pipeline.log                  # Processing log
```

## Method Details

### Reference-Based Approach (for alignments >100 sequences)
- Uses actual reference sequence from Step 02 when available
- Falls back to consensus sequence if reference not found
- O(n·L) complexity instead of O(n²) for pairwise
- Proper site normalization using syn_fraction_by_pos
- No per-sequence normalization to avoid bias

### Pairwise Approach (for alignments ≤100 sequences)
- Full pairwise comparisons with NG86 method
- Aggregates rates (sum of dN/dS) instead of averaging omegas
- More accurate for small sample sizes

### Site Normalization
- Calculates synonymous opportunity at each codon position
- Accounts for codon degeneracy and genetic code structure
- Cached computation for efficiency (LRU cache with 4096 entries)

## Key Metrics

### dN/dS Ratio Interpretation
- **< 0.1**: Strong purifying selection
- **0.1-0.5**: Purifying selection
- **0.5-1.0**: Weak purifying selection
- **≈ 1.0**: Neutral evolution
- **> 1.0**: Positive selection

### Site-Based SNP Metrics (NEW)
- **codon_variable_sites**: Total variable codon positions
- **codon_syn_sites**: Variable positions where all variants encode same amino acid
- **codon_nonsyn_sites**: Variable positions that change amino acid
- **codon_syn_fraction**: Proportion of variable sites that are synonymous

### Additional Metrics
- **S_sites**: Number of synonymous sites (opportunity)
- **N_sites**: Number of non-synonymous sites (opportunity)
- **ng86_total_codons**: Total valid codons analyzed (fixed calculation)
- **Variable sites**: Positions with sequence variation
- **Parsimony informative**: Sites useful for phylogeny
- **Gap percentage**: Alignment quality indicator

## Statistical Analysis

The pipeline performs:
- **Mann-Whitney U test**: Non-parametric comparison of dN/dS distributions
- **Effect size**: Rank-biserial correlation
- **Descriptive statistics**: Mean, median, quartiles
- **Selection pressure classification**: Categorization by dN/dS ranges

## Resource Requirements

- **Time**: 1-2 hours (depending on dataset size and mode)
- **Memory**: 32GB
- **CPUs**: 20 (adjustable with --threads)
- **Partition**: cpu

## Dependencies

Required Python packages:
- BioPython (≥1.79)
- pandas (≥1.3)
- numpy (≥1.21)
- matplotlib (≥3.4)
- seaborn (≥0.11)
- scipy (≥1.7)

## Examples

### Analyze only operon genes with 10 threads
```bash
python dnds_pipeline.py \
    --mode operon \
    --operon-dir ../05_operon_assembly_extraction/output/msa/dna_alignments \
    --output-dir operon_only_output \
    --threads 10
```

### Compare methods on core genes
```bash
python dnds_pipeline.py \
    --mode core \
    --core-dir ../04_core_gene_analysis/output/core_gene_alignments \
    --output-dir core_comparison \
    --include-simple \
    --threads 20
```

### Submit SLURM job for operon analysis only
```bash
sbatch run_dnds_analysis.sh --mode operon
```

## Troubleshooting

If the pipeline fails:
1. Check input directories exist and contain alignments
2. Review the pipeline log: `output/dnds_pipeline.log`
3. Check SLURM error logs: `dnds_pipeline_*.err`
4. Verify environment activation and package installation
5. Ensure sufficient memory for large datasets
6. For timeout issues, consider using `--mode` to analyze subsets

## Performance Tips

- Use `--mode operon` or `--mode core` to analyze subsets faster
- Alignments >100 sequences automatically use the faster reference-based method
- The `--include-simple` flag adds significant computation time
- Increase `--threads` for faster parallel processing
- Cache is automatically managed (4096 codon entries)