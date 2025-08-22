# Step 6: Comparative Diversity Analysis

This directory contains scripts for comparative analysis of conservation patterns between core genes and operon genes.

## Overview

This step creates conservation ranking visualizations comparing operon genes to core genes using pre-computed metrics from steps 04 and 05. The analysis focuses on visual comparisons showing where operon genes rank among all core genes.

## Input Data

### Exact Files Used

The analysis reads these specific pre-computed CSV files:

1. **Core Gene Metrics** (from Step 04):
   - File: `../04_core_gene_analysis/output/core_gene_conservation_metrics.csv`
   - Contains: Conservation scores and pairwise identity for 1,280 core genes
   - Columns used: `gene`, `mean_conservation`, `mean_pairwise_identity`

2. **Operon Gene Metrics** (from Step 05):
   - Primary file: `../05_operon_assembly_extraction/output/mappings/aa_nt_mapping/assemblies/msa/operon_conservation_metrics.csv`
   - Alternative (if primary not found): `../05_operon_assembly_extraction/output/mappings/aa_nt_mapping/prokka/msa/operon_conservation_metrics.csv`
   - Contains: Conservation scores for 7 operon genes (frpC, glpC, ptsA, ptsB, ptsC, ptsD, fruR)
   - Columns used: `gene`, `mean_conservation`, `mean_pairwise_identity`

**Note**: The script uses Strategy D (assembly-based tblastn) metrics as the primary source because it provides the most authentic sequence data from raw assemblies.

## Output Files

All outputs are written to `output/`:
- `conservation_ranking_bars.png/pdf` - Three-panel bar plot showing top 20, operon genes, and bottom 20 conserved genes
- `conservation_distribution.png/pdf` - Histogram of conservation distribution with operon genes highlighted
- `operon_conservation_summary.csv` - Summary table with operon gene rankings and percentiles

## Scripts

### `comparative_analysis.py`
Main comparative analysis script that:
- Loads pre-computed conservation metrics from steps 04 and 05
- Creates conservation ranking bar plots
- Calculates operon gene rankings among core genes
- Generates summary statistics
- Does NOT re-calculate MSAs or conservation scores

### `run_analysis.sh`
SLURM script to run the comparative analysis pipeline.

## Usage

```bash
# Run comparative analysis via SLURM
sbatch run_analysis.sh

# Or run directly
python comparative_analysis.py
```

## Visualization Details

### Conservation Ranking Bar Plot
Three-panel visualization showing:
1. **Top 20 Panel**: Most conserved core genes (green bars)
2. **Operon Panel**: All 7 operon genes with their rankings (red bars)
3. **Bottom 20 Panel**: Least conserved core genes (dark red bars)

Each bar shows:
- Conservation percentage
- Rank among all core genes
- Percentile position

### Conservation Distribution
Histogram showing the distribution of conservation scores across all core genes with operon genes highlighted as vertical lines.

## Key Outputs

The analysis generates:
- Visual ranking of operon genes among core genes
- Percentile positions showing relative conservation
- Summary statistics comparing mean conservation levels
- Individual gene rankings and conservation scores

## Dependencies

This step requires completed analyses from:
- Step 04: Core gene analysis must be complete
- Step 05: Operon extraction and MSA must be complete

## Notes

- Analysis uses pre-computed metrics to avoid redundant calculations
- Statistical tests are non-parametric to handle non-normal distributions
- Multiple extraction strategies from step 05 are compared simultaneously
- Results provide quantitative support for manuscript conclusions