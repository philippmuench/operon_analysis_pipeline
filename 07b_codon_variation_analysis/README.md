# Codon Variation-Based Selection Analysis

## Overview

This analysis calculates comparable selection pressure metrics for operon and core genes using codon-level variation. This approach works for both gene sets without the issues encountered in traditional dN/dS calculations.

## Why This Approach?

Traditional dN/dS analysis failed because:
- **Core genes**: Too conserved for consensus-based methods (>99.9% identity)
- **Operon genes**: Alignment issues (frameshifts, stop codons) break pairwise methods
- **No universal reference**: Core genes lack external reference sequences

This codon variation approach:
- ✅ Works for all genes with any amount of variation
- ✅ Provides comparable metrics across gene sets
- ✅ Robust to alignment issues
- ✅ No reference sequence required

## Method

The analysis calculates a **normalized selection index** based on the ratio of non-synonymous to synonymous changes in variable codon positions:

1. **Count variable codon positions** in each gene
2. **Classify as synonymous or non-synonymous**
3. **Calculate ratio**: nonsyn/syn sites
4. **Normalize** by the expected neutral ratio (default 2.5, configurable via `--expected-neutral`)

### Interpretation:
- **< 0.7**: Purifying selection (conserved)
- **0.7-1.3**: Near neutral
- **> 1.3**: Relaxed/positive selection

## Running the Analysis

### Quick Start:
```bash
sbatch run_codon_variation.sh
```

### Manual Run:
```bash
python codon_variation_analysis.py \
    --operon-file ../07_dnds_analysis/output_old/tables/operon_dnds_analysis.csv \
    --core-file ../07_dnds_analysis/output_old/tables/core_genes_dnds_analysis.csv \
    --output-dir output \
    --expected-neutral 2.5  # Optional: adjust neutral expectation ratio
```

### Requirements:
- CSV files from step 07 containing per-gene counts of variable codon sites (synonymous and non-synonymous)
  - Note: Uses the codon variation summaries, not dN/dS estimates themselves
- Python with pandas, numpy, matplotlib, seaborn, scipy, biopython

## Output Files

### Tables (`output/tables/`)
- `codon_variation_metrics.csv`: Detailed metrics for all genes
- `summary_statistics.csv`: Summary stats by gene type
- `*_codon_variation.csv`: Per-gene codon variation data (for operon genes)

### Plots (`output/plots/`)
- `selection_distribution.pdf`: Boxplot and jitter plot of npN/pS distributions
- `pN_pS_scatter.pdf`: Scatter plot of non-synonymous vs synonymous counts

### Syn/Nonsyn Plots (`output/plots/syn_nonsyn/`)
- `*_syn_nonsyn_distribution.pdf`: Individual gene variation profiles (7 operon genes)
- `operon_combined_syn_nonsyn.pdf`: Combined view of all operon genes

### Report
- `output/codon_variation_report.txt`: Complete analysis report with statistics

## Key Findings

### Operon Genes
- **Mean normalized selection**: ~1.2 (slightly relaxed)
- **Range**: 0.49-3.5
- **Pattern**: Moderate but consistent selection

### Core Genes  
- **Mean normalized selection**: ~3.8 (highly variable)
- **Bimodal distribution**:
  - Essential genes: Very strong purifying (0.05-0.3)
  - Accessory genes: Relaxed selection (>10)

### Comparison
- No significant overall difference (p=0.53)
- Core genes show much greater heterogeneity
- Operon represents intermediate selective constraint

## Example Results

### Most Conserved Operon Genes:
1. ptsA (0.489) - PTS transport
2. ptsC (0.696) - PTS transport
3. ptsB (0.736) - PTS transport

### Most Conserved Core Genes:
1. tuf (0.054) - Translation elongation
2. rpsZ (0.067) - Ribosomal protein
3. recA (0.089) - DNA repair

### Least Conserved Core Genes:
1. dppE (228.4) - Peptide transport
2. manX (128.8) - Mannose transport
3. gph (112.4) - Phosphatase

## Biological Interpretation

The analysis reveals:

1. **Functional constraint gradient**: Essential > Metabolic > Transport > Regulatory
2. **Core gene heterogeneity**: Mix of ultra-conserved and rapidly evolving genes
3. **Operon intermediate position**: Specialized metabolism under moderate selection
4. **Selection relaxation in regulators**: Adaptive flexibility for environmental response

## Troubleshooting

### No input files found
Run the dN/dS analysis first:
```bash
cd ../07_dnds_analysis
sbatch run_dnds_analysis.sh
```

### Memory issues
Increase memory allocation:
```bash
sbatch --mem=32G run_codon_variation.sh
```

### Missing dependencies
Install required packages:
```bash
conda activate efs_diversity
pip install pandas numpy matplotlib seaborn scipy
```