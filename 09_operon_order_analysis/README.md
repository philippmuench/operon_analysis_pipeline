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
- Directory: `../03_blast_search/output/aa_nt_mapping/assemblies/`
- Contains genomic positions of operon genes from BLAST hits

## Canonical Operon Order

Reference gene order: `frpC → glpC → ptsD → ptsC → ptsB → ptsA → fruR`

## Scripts

- `analyze_operon_order.py` - Analyzes gene positions and order
- `visualize_operon_synteny.py` - Creates visualization plots
- `manuscript_numbers.py` - Generates statistics for manuscript
- `run_operon_order_analysis.sh` - SLURM submission script (includes all steps)
- `run_manuscript_stats.sh` - Generate manuscript statistics only

## Usage

```bash
# Run complete analysis (including manuscript stats)
sbatch run_operon_order_analysis.sh

# Generate manuscript statistics only (if analysis already done)
sbatch run_manuscript_stats.sh
```

## Output Files

Located in `output/`:
- `operon_order_analysis.csv` - Gene order and organization for each genome
- `operon_order_summary.txt` - Summary statistics
- `manuscript_stats.txt` - Key statistics formatted for manuscript
- `operon_synteny_examples.png` - Visual examples of different organizations
- `gene_presence_heatmap.png` - Gene presence/absence patterns
- `operon_organization_summary.png` - Statistical summary plots

## Operon Types

Genomes are classified as:
- **Complete and conserved**: All 7 genes in canonical order
- **Complete with rearrangement**: All 7 genes but reordered
- **Complete but fragmented**: All 7 genes across multiple contigs
- **Near-complete**: 5-6 genes present
- **Partial**: 3-4 genes present
- **Minimal**: 1-2 genes present