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

Uses multiple data sources:
- **Prokka annotations** from Step 01: `../01_prokka_annotation/output/prokka_results/`
- **Reference genes** from Step 02: `../02_reference_operon_extraction/output/operon_genes_nt.fasta`
- **BLAST mappings** from Step 03: `../03_blast_search/output/blast_results_prokka_variants/`
  - Contains gene-to-gene mappings (reference operon → Prokka predicted genes)
  - Essential for linking Prokka locus tags to operon gene names

## Scripts

- `analyze_start_sites.py` - Main analysis of start codons and RBS
  - Parses Prokka GFF/FAA files for gene sequences
  - Uses BLAST results to map Prokka locus tags to operon gene names
  - Detects RBS motifs and alternative start sites
  - Outputs detailed TSV with results per gene instance
- `manuscript_numbers.py` - Generates statistics for manuscript
  - Processes the TSV output from analyze_start_sites.py
  - Calculates summary statistics per gene and overall
  - Formats results for manuscript text
- `run_start_site_analysis.sh` - SLURM submission script
  - Runs complete pipeline: analysis + manuscript stats
  - Automatically activates conda environment
  - Handles all path configurations

## Usage

```bash
# Run complete analysis (including manuscript stats)
cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/08_start_site_analysis
sbatch run_start_site_analysis.sh

# Check job status
squeue -u $USER

# Monitor output
tail -f start_site_*.out
```

### Prerequisites

Before running, ensure previous steps are complete:
1. Step 01 (Prokka annotation) must be complete
2. Step 02 (Reference extraction) must be complete  
3. Step 03 (BLAST search) must have `blast_results_prokka_variants/` populated

### Expected Runtime

- ~24 hours for full dataset (8,587 genomes)
- Processes genomes in parallel using 8 CPUs
- Memory usage: ~8GB

## Output Files

Located in `output/`:
- `start_site_summary.tsv` - Detailed results for each gene instance
- `manuscript_stats.txt` - Key statistics formatted for manuscript
- `report/report.html` - Visual report with highlighted start sites and RBS

## Key Metrics

- **Start codon distribution**: % using ATG vs GTG vs TTG
  - Most genes use ATG (canonical)
  - ptsA uniquely uses TTG (87% of instances)
- **Alternative starts**: Genes with potential upstream start sites
  - ~87% of genes have upstream ATG codons
  - ~72% have upstream ATG with associated RBS
- **RBS presence**: Genes with detected Shine-Dalgarno sequences
  - AGGAGG motif and variants (GGAGG, AGGAG, GAGG)
- **RBS spacing**: Distance between RBS and start codon
  - Optimal: 5-9 nt
  - Acceptable: 4-13 nt

## Classifications

Genes are classified based on their start sites:
- **Alt_upstream_with_RBS**: Alternative upstream start with RBS present
  - Suggests possible misannotation or translational flexibility
- **Upstream_no_RBS**: Upstream start exists but lacks RBS
  - Less likely to be functional
- **No_upstream**: No alternative upstream start found
  - Annotation likely correct
- **Partial_gene**: Gene appears truncated
  - Missing expected start region
- **Same_start**: Predicted start matches expectation
  - High confidence in annotation

## Biological Significance

The analysis reveals:
- **ptsA's unique TTG usage**: Suggests specialized translational regulation
- **High alternative start prevalence**: May indicate annotation challenges or biological flexibility
- **Operon-wide patterns**: Consistent start codon usage within functional modules