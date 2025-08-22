# Step 9: Operon Order and Synteny Analysis

Analyzes how operon genes are organized and ordered in bacterial genomes.

## Overview

This step examines:
- Which operon genes are present in each genome
- The order and orientation of genes
- Whether genes are on the same contig and strand
- How well the gene order is conserved

## Input Data

Uses BLAST results from Step 03:
- Directory: `../03_blast_search/output/blast_results/`
- Contains genomic positions of operon genes from BLAST hits

## Canonical Operon Order

Reference gene order: `frpC → glpC → ptsD → ptsC → ptsB → ptsA → fruR`

## Scripts

- `operon_order_pipeline.py` - Unified pipeline script with all analysis functionality
- `run_pipeline.sh` - SLURM submission script to run the complete pipeline

## Usage

```bash
# Run complete pipeline (analysis, visualization, and statistics)
sbatch run_pipeline.sh

# Or run specific steps interactively
python operon_order_pipeline.py --steps analyze          # Run analysis only
python operon_order_pipeline.py --steps visualize        # Create visualizations only
python operon_order_pipeline.py --steps stats            # Generate statistics only
python operon_order_pipeline.py --steps all              # Run all steps (default)

# With custom parameters
python operon_order_pipeline.py \
    --blast-dir ../03_blast_search/output/blast_results \
    --output-dir output \
    --min-identity 80 \
    --steps all
```

## Output Files

Located in `output/`:
- `operon_order_analysis.csv` - Detailed gene order and organization for each genome
- `operon_order_summary.txt` - Summary statistics
- `manuscript_stats.txt` - Comprehensive statistics formatted for manuscript
- `gene_cooccurrence_matrix.csv` - Gene co-occurrence matrix (CSV only, no heatmap)
- `non_canonical_gene_orders.txt` - List of genomes with non-canonical gene orders
- `operon_organization_types.png` - **Horizontal bar chart** of operon type distribution
- `gene_order_visualization.pdf` - Visual examples showing:
  - Canonical reference order
  - Different organization types with genome examples
  - Fragmentation indicators (`//` for genes on different contigs)
  - Missing gene counts for partial operons

## Operon Types

Genomes are classified as:
- **Complete and conserved**: All 7 genes in canonical order on same contig
  - "Conserved" specifically means maintaining the canonical gene order
- **Complete with rearrangement**: All 7 genes but in different order
- **Complete but fragmented**: All 7 genes across multiple contigs
- **Near-complete**: 5-6 genes present
- **Partial**: 3-4 genes present
- **Minimal**: 1-2 genes present
- **No operon**: 0 genes detected

## Visualization Improvements

### Gene Co-occurrence Heatmap Removed
The heatmap was removed because all values were ~1.0 (uninformative). When the operon 
is present, all genes co-occur nearly 100% of the time. The CSV matrix is still 
generated for reference.

### Enhanced Visualizations
- **Bar chart instead of pie chart**: Eliminates overlapping labels
- **Fragmentation indicators**: `//` symbol shows where genes span contigs
- **Missing gene counts**: Shows "(X missing)" for partial operons
- **Consistent gene sizing**: All gene boxes use same width for easier comparison