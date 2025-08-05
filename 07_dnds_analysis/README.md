# Step 6: dN/dS and Substitution Analysis

This directory contains scripts for analyzing selection pressure and substitution patterns in the operon genes.

## Overview

This step performs specialized analysis of:
- **dN/dS ratios** (ratio of non-synonymous to synonymous substitutions)
- **Selection pressure** (purifying vs positive selection)
- **Substitution patterns** (transitions vs transversions)
- **Codon usage** and evolutionary constraints

## Input Data

Uses data from previous pipeline steps:
- Multiple sequence alignments from `../05_diversity_analysis/output/msa/`
- BLAST results from `../03_blast_search/output/`

## Main Scripts

### **dN/dS Analysis Scripts**
- `analyze_dnds.py` - Calculate dN/dS ratios for gene alignments
- `compare_dnds_results.py` - Compare dN/dS across different genes
- `run_dnds_analysis.sh` - Main dN/dS analysis pipeline

### **Substitution Analysis Scripts**
- `analyze_operon_substitutions.py` - Analyze substitution patterns in operons
- `analyze_substitution_types.py` - Classify substitution types (Ti/Tv)
- `compare_substitutions_all.py` - Comprehensive substitution comparison
- `submit_substitution_comparison.sh` - Submit substitution analysis job

## Usage

### **Run dN/dS Analysis**
```bash
# Submit dN/dS analysis job
sbatch run_dnds_analysis.sh
```

### **Run Substitution Analysis**
```bash
# Submit substitution comparison job
sbatch submit_substitution_comparison.sh
```

## Output Files

Results are typically saved to an `output/` directory within this folder:

### dN/dS Analysis Results
- `dnds_results.csv` - dN/dS ratios for each gene
- `selection_summary.txt` - Selection pressure classification
- `dnds_comparison_plots.png` - Comparative visualizations

### Substitution Analysis Results
- `substitution_patterns.csv` - Detailed substitution classifications
- `ti_tv_ratios.csv` - Transition/transversion ratios
- `substitution_summary_plots.png` - Pattern visualizations

## Key Metrics Explained

### dN/dS Ratio
- **dN/dS < 0.5**: Strong purifying selection (negative selection)
- **dN/dS ≈ 1.0**: Neutral evolution
- **dN/dS > 1.0**: Positive selection

### Ti/Tv Ratio
- **Expected ~2.0**: For neutral evolution
- **Higher values**: Indicate selection constraints
- **Lower values**: May indicate relaxed selection

## Prerequisites

Before running these analyses, ensure you have completed:
1. Step 3: BLAST search
2. Step 5: Diversity analysis with MSA creation

## Dependencies

- Python 3.x with BioPython
- Multiple sequence alignments from Step 5
- BLAST results from Step 3