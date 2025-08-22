# Step 8: Start Site Analysis

Analyzes translation initiation sites and ribosome binding sites for operon genes.

## Overview

This step examines:
- Start codon usage (ATG, GTG, TTG)
- Alternative upstream start sites
- Ribosome binding sites (Shine-Dalgarno sequences)
- Discrepancies between predicted and potential start sites
- **NEW: Metadata stratification** by source niche and country
- **NEW: Laboratory adaptation patterns** in start codon usage

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
  - **NEW: Stratifies results by metadata** (source niche, country)
  - **NEW: Creates gene-specific visualizations** for all 7 operon genes
  - Outputs detailed TSV with results per gene instance
- `manuscript_numbers.py` - Generates statistics for manuscript
  - Processes the TSV output from analyze_start_sites.py
  - Calculates summary statistics per gene and overall
  - **NEW: Statistical analysis of ptsA laboratory pattern** (Fisher's exact test)
  - Formats results for manuscript text
- `run_start_site_analysis.sh` - SLURM submission script
  - Runs complete pipeline: analysis + manuscript stats
  - **NEW: Supports visualization-only mode** for quick plot regeneration
  - Automatically activates conda environment
  - Handles all path configurations

## Usage

```bash
# Run complete analysis (including manuscript stats)
cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/08_start_site_analysis
sbatch run_start_site_analysis.sh

# Run visualization-only mode (uses existing results)
sbatch run_start_site_analysis.sh visualize

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
- **NEW: Stratification files**:
  - `start_codon_by_source_niche.tsv` - Start codon usage by source niche
  - `start_codon_by_country.tsv` - Start codon usage by country
  - `{gene}_by_source_niche.tsv` - Gene-specific niche analysis (7 files)
  - `{gene}_by_country.tsv` - Gene-specific country analysis (7 files)
  - `metadata_stratification_summary.txt` - Summary of stratified results
- **NEW: Visualization plots** in `plots/`:
  - `{gene}_by_country.png/pdf` - Stacked bar charts for each gene by country
  - `{gene}_by_source_niche.png/pdf` - Stacked bar charts for each gene by niche
  - `start_codon_by_country_detailed.png/pdf` - Overall country comparison
  - `ttg_usage_heatmap.png/pdf` - TTG usage patterns across genes and niches

## Key Metrics

- **Start codon distribution**: % using ATG vs GTG vs TTG
  - Most genes use ATG (87.3% overall)
  - ptsA uniquely uses TTG (87.3% of instances)
  - GTG essentially absent (0.03%)
- **Alternative starts**: Genes with potential upstream start sites
  - 86.9% of genes have upstream start codons
  - 72.3% have upstream starts with associated RBS
- **RBS presence**: Genes with detected Shine-Dalgarno sequences
  - AGGAGG motif and variants (GGAGG, AGGAG, GAGG)
- **RBS spacing**: Distance between RBS and start codon
  - Mean: 30 nt, Median: 24 nt
  - Optimal: 5-9 nt, Acceptable: 4-13 nt
- **Laboratory adaptation** (NEW):
  - Laboratory strains: 55.1% ATG usage in ptsA
  - Other niches: 10.9% ATG usage in ptsA
  - Odds ratio: 10.0 (p < 0.001, Fisher's exact test)

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
- **Laboratory adaptation pattern**: Laboratory strains show 10-fold higher ATG usage in ptsA
  - Indicates adaptation to cultivation conditions affects translation initiation
  - Statistical significance: p = 4.4×10⁻⁷⁹ (Fisher's exact test)
  - May reflect selection for altered PTS component stoichiometry
- **High alternative start prevalence**: May indicate annotation challenges or biological flexibility
  - 86.9% of genes have potential upstream starts
  - Most (72.3%) have associated RBS motifs
- **Operon-wide patterns**: Consistent start codon usage within functional modules
  - Six genes exclusively use ATG (99-100%)
  - Only ptsA shows TTG preference

## Metadata Stratification

Analysis includes stratification by:
- **Source Niche**: Human, Animal, Environmental, Food, Laboratory
- **Country**: Geographic origin of isolates
- **Source Details**: Specific isolation sources

Key finding: Laboratory strains show distinct adaptation in ptsA start codon usage,
suggesting cultivation conditions select for altered translational regulation.