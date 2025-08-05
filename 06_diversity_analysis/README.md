# Step 6: Diversity Analysis

This directory contains scripts for analyzing sequence diversity and conservation patterns in both core genes and operon genes.

## Overview

This step performs comprehensive diversity analysis using MSAs from both core genes (step 04) and operon genes (step 05) to compare conservation patterns and evolutionary pressures.

## Input Data

- Core gene MSAs from `../04_core_gene_analysis/output/core_gene_alignments/`
- Operon gene MSAs from `../05_operon_extraction/output/msa/`
- BLAST results from `../03_blast_search/output/`

## Output Files

All outputs are written to the `output/` directory:
- `diversity_results.csv`: Comprehensive diversity metrics
- `conservation_profiles.png`: Position-based conservation plots
- `promoter_conservation_with_pribnow.png`: Promoter analysis with Pribnow box
- `gap_analysis.png`: Gap pattern analysis
- `diversity_summary.png`: Overall summary plots

## Scripts

### Main Analysis Scripts

#### `analyze_diversity.py`
Main diversity analysis script for coding genes.
- Calculates pairwise sequence identity
- Shannon entropy-based conservation scores
- Position-wise conservation analysis
- Creates conservation profile plots

#### `analyze_enhanced_diversity.py`
Enhanced analysis for both coding and non-coding sequences.
- Comprehensive metrics calculation
- Multiple visualization types
- Statistical summaries

#### `analyze_msa_diversity.py`
Specialized MSA-based diversity analysis.
- Works directly with alignment files
- Advanced diversity metrics
- Phylogenetic considerations

### Comparison Scripts

#### `analyze_all_core_genes.py`
Analysis focusing on core gene patterns.
- Core gene specific metrics
- Comparison with operon genes
- Baseline establishment

#### `calculate_core_gene_diversity.py`
Dedicated core gene diversity calculator.
- Statistical analysis of core genes
- Distribution analysis
- Outlier detection

### Visualization Scripts

#### `generate_diversity_summary.py`
Creates comprehensive summary visualizations.
- Comparative plots
- Statistical distributions
- Publication-ready figures

### Pipeline Scripts

#### `run_complete_diversity_analysis.sh`
Master pipeline script that runs complete diversity analysis.
- Orchestrates all analysis steps
- Handles both core and operon genes
- Creates all plots and summaries

#### `run_diversity_analysis_array.sh`
Array job version for large-scale analysis.

## Usage

### Complete Analysis
```bash
# Run full diversity analysis pipeline
sbatch run_complete_diversity_analysis.sh
```

### Individual Analyses
```bash
# Analyze coding genes
python analyze_diversity.py

# Analyze with enhanced features
python analyze_enhanced_diversity.py

# Focus on core genes
python analyze_all_core_genes.py
```

### Custom Analysis
```bash
# Specific MSA analysis
python analyze_msa_diversity.py --input msa_directory/ --output results/
```

## Key Metrics

### Conservation Metrics
- **Shannon entropy**: Position-wise conservation scores
- **Pairwise identity**: Sequence similarity measurements
- **Gap analysis**: Insertion/deletion patterns
- **Conservation windows**: Regional conservation patterns

### Diversity Metrics
- **Nucleotide diversity (π)**: Average pairwise differences
- **Segregating sites**: Variable positions
- **Transition/transversion ratios**: Mutation pattern analysis
- **Selection signatures**: Evidence for natural selection

### Comparison Metrics
- **Core vs. operon**: Comparative conservation analysis
- **Gene-specific patterns**: Individual gene diversity profiles
- **Functional categories**: Analysis by gene function

## Key Features

- **Dual comparison**: Core genes vs. operon genes
- **Multiple metrics**: Comprehensive diversity measurement
- **Visualization**: Publication-ready plots
- **Statistical analysis**: Significance testing
- **Scalable**: Handles large gene sets efficiently

## Integration

This step uses output from:
- `../04_core_gene_analysis/`: Core gene MSAs and metrics
- `../05_operon_extraction/`: Operon gene MSAs

This step provides input for:
- `../07_dnds_analysis/`: dN/dS ratio calculations
- Publication figures and tables

## Scientific Rationale

The comparison between core genes and operon genes allows:
1. **Baseline establishment**: Core genes provide "normal" conservation levels
2. **Pathway analysis**: Operon genes reveal pathway-specific evolution
3. **Selection analysis**: Different evolutionary pressures identification
4. **Functional insights**: Conservation-function relationships

## Notes

- Analysis requires high-quality MSAs from previous steps
- Statistical significance depends on number of sequences analyzed
- Core gene analysis provides crucial baseline for operon comparison
- Results directly inform biological interpretation of operon evolution