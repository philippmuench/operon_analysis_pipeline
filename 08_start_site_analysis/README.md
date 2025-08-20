# Step 8: Start Site Analysis

Analyzes translation initiation sites and ribosome binding sites for operon genes.

## Overview

This step examines:
- Start codon usage (ATG, GTG, TTG)
- Alternative upstream start sites
- Ribosome binding sites (Shine-Dalgarno sequences)
- Discrepancies between predicted and potential start sites

## Unique Analysis

This step provides analysis NOT covered in other steps:
- **Translation initiation**: Start codon preferences
- **RBS detection**: Shine-Dalgarno motif identification
- **Gene boundary validation**: Alternative start site detection
- **Regulatory signals**: RBS spacer distance analysis

## Input Data

Uses Prokka annotations from Step 01:
- Directory: `../01_prokka_annotation/output/prokka_results/`
- Reference genes: `../02_reference_operon_extraction/output/operon_genes_nt.fasta`

## Scripts

- `analyze_start_sites.py` - Main analysis of start codons and RBS
- `manuscript_numbers.py` - Generates statistics for manuscript
- `run_start_site_analysis.sh` - SLURM submission script (includes all steps)
- `run_manuscript_stats.sh` - Generate manuscript statistics only

## Usage

```bash
# Run complete analysis (including manuscript stats)
sbatch run_start_site_analysis.sh

# Generate manuscript statistics only (if analysis already done)
sbatch run_manuscript_stats.sh
```

## Output Files

Located in `output/`:
- `start_site_summary.tsv` - Detailed results for each gene instance
- `manuscript_stats.txt` - Key statistics formatted for manuscript
- `report/report.html` - Visual report with highlighted start sites and RBS

## Key Metrics

- **Start codon distribution**: % using ATG vs GTG vs TTG
- **Alternative starts**: Genes with potential upstream start sites
- **RBS presence**: Genes with detected Shine-Dalgarno sequences
- **RBS spacing**: Distance between RBS and start codon (typically 4-13 nt)

## Classifications

Genes are classified based on their start sites:
- **Alt_upstream_with_RBS**: Alternative upstream start with RBS present
- **Upstream_no_RBS**: Upstream start exists but lacks RBS
- **No_upstream**: No alternative upstream start found
- **Partial_gene**: Gene appears truncated
- **Same_start**: Predicted start matches expectation