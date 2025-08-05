# Step 5: Complete Diversity Analysis

This directory contains scripts for comprehensive sequence conservation and diversity analysis of operon genes and non-coding regions.

## Overview

This step performs complete analysis of sequence diversity across the 48 genomes that contain the complete operon:

### Coding Gene Analysis:
- Pairwise sequence identity
- Conservation scores per position
- Substitution patterns (transitions vs transversions)
- Selection pressure (dN/dS ratios)
- Gap analysis in alignments

### Non-coding Region Analysis:
- Promoter sequence conservation
- Promoter identity across genomes
- Comparative conservation plots

## Input Data

Uses data from previous pipeline steps:
- BLAST results from `../03_blast_search/output/`
- Prokka annotations from `../01_prokka_annotation/output/prokka_results/`

## Main Script

### **Complete Pipeline**
- `run_complete_diversity_analysis.sh` - **ONE SCRIPT THAT RUNS EVERYTHING**

## Usage

### **Run Complete Analysis**
```bash
# Submit the complete analysis job
sbatch run_complete_diversity_analysis.sh
```

## Output Files

All outputs are saved to the `output/` directory:

### Sequences
- `output/operon_sequences/` - FASTA files for each operon gene

### Alignments
- `output/msa/dna_alignments/` - DNA multiple sequence alignments
- `output/msa/protein_alignments/` - Protein sequences and alignments
- `output/msa/alignment_summary.csv` - Summary of alignment statistics

### Analysis Results
- `output/diversity_analysis/diversity_results.csv` - Detailed diversity metrics
- `output/diversity_analysis/conservation_profiles.png` - Conservation plots for each gene
- `output/diversity_analysis/diversity_summary.png` - Summary visualization
- `output/diversity_analysis/diversity_summary.txt` - Text summary report

## Key Findings (from test run)

- **Overall conservation**: >99% DNA identity, >99% protein identity
- **Most conserved gene**: ptsA (99.57% DNA, 99.92% protein)
- **Least conserved gene**: glpC (98.92% DNA, 98.73% protein)
- **Selection**: Strong purifying selection on all genes (dN/dS < 0.5)
- **All operon genes are highly conserved**, indicating strong functional constraints

## Metrics Explained

### Conservation Score
- 1.0 = Completely conserved position
- 0.0 = Highly variable position
- Based on Shannon entropy

### dN/dS Ratio
- dN/dS < 0.5: Strong purifying selection
- dN/dS = 1.0: Neutral evolution
- dN/dS > 1.0: Positive selection

### Ti/Tv Ratio
- Expected ~2.0 for neutral evolution
- Higher values indicate selection constraint

## Individual Python Scripts (Advanced Users)

The main pipeline uses these individual scripts in sequence:

### Core Analysis Scripts
- `extract_operon_sequences.py` - Extracts coding gene sequences
- `extract_noncoding_sequences.py` - Extracts promoter sequences  
- `create_msa.py` - Creates MSAs for coding sequences
- `create_enhanced_msa.py` - Creates MSAs for non-coding sequences
- `analyze_diversity.py` - Analyzes coding gene diversity
- `analyze_enhanced_diversity.py` - Analyzes all sequence types
- `create_gap_analysis.py` - Creates gap analysis plots
- `create_promoter_msa_from_blast.py` - Creates promoter plots from BLAST data

### Additional Analysis Scripts
- `calculate_core_gene_diversity.py` - Core gene analysis
- `generate_diversity_summary.py` - Summary generation
- `analyze_msa_diversity.py` - MSA-specific analysis
- Various other specialized diversity analysis scripts

**Note**: dN/dS and substitution analysis scripts have been moved to `../06_dnds_analysis/`

## What the Complete Pipeline Does

### Step 1: Extract Sequences
- Coding genes from BLAST results
- Promoter sequences from non-coding BLAST results

### Step 2: Create Alignments  
- Multiple sequence alignments for all genes
- Promoter sequence alignments

### Step 3: Analyze Diversity
- Conservation scores and identity metrics
- dN/dS ratios and selection pressure
- Gap patterns in alignments

### Step 4: Generate Visualizations
- Conservation profiles for each gene
- Promoter conservation plots
- Gap frequency plots
- Comparative summary plots

### Promoter Data Location
The `has_promoter` information is available in:
```
../03_blast_search/output/operon_simple_summary.csv
```

Columns of interest:
- `has_promoter`: 1 if promoter found, 0 otherwise
- `promoter_identity`: BLAST identity percentage for promoter
- `promoter_coverage`: BLAST coverage percentage for promoter