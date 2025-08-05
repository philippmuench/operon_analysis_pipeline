# Step 5: Diversity Analysis

This directory contains scripts for analyzing sequence conservation and selection pressure.

## Analysis Components

### 1. Core Gene Diversity
Calculate pairwise identity for all core genes to establish baseline conservation levels.

```bash
python calculate_core_gene_diversity.py \
    --input_dir ../../core_gene_sequences_95pct \
    --output ../../results/core_gene_diversity.csv
```

### 2. MSA-based Diversity
Analyze multiple sequence alignments for detailed conservation metrics.

```bash
python analyze_msa_diversity.py \
    --msa_dir ../../msa_output \
    --output ../../results/msa_diversity.csv
```

### 3. dN/dS Analysis
Calculate synonymous vs non-synonymous substitution rates to assess selection pressure.

```bash
python analyze_dnds.py \
    --input_dir ../../msa_output \
    --output ../../results/dnds_results.csv
```

### 4. Comparative Analysis
Compare operon genes against core genome background.

```bash
python compare_substitutions_all.py \
    --operon_dir ../../msa_output \
    --core_dir ../../core_gene_sequences_95pct \
    --output ../../results/comparison.csv
```

## Run Complete Pipeline

```bash
sbatch run_diversity_pipeline.sh
```

This will:
1. Calculate diversity for all core genes
2. Analyze operon gene conservation
3. Compute dN/dS ratios
4. Compare operon vs core genes
5. Generate summary statistics

## Output Files

- `core_gene_diversity.csv` - Conservation metrics for each core gene
- `operon_diversity_detailed.csv` - Detailed operon conservation
- `dnds_results.csv` - dN/dS ratios and selection analysis
- `diversity_summary.json` - Comprehensive summary statistics

## Key Metrics

### Conservation (Pairwise Identity)
- **>98%**: Highly conserved
- **95-98%**: Conserved
- **90-95%**: Moderately conserved
- **<90%**: Variable

### Selection (dN/dS Ratio)
- **<0.5**: Strong purifying selection
- **0.5-1.0**: Weak purifying selection
- **~1.0**: Neutral evolution
- **>1.5**: Positive selection

## Results Summary

- Operon genes: 98.7-99.6% identity (highly conserved)
- Core genes: 94.7% ± 9.1% average identity
- Operon dN/dS: <0.5 (strong purifying selection)
- All operon genes rank in top 30% for conservation