# Analysis Re-run Plan

## Overview

Re-running the operon analysis pipeline to include 23 new genomes from `000_new_genomes/`. The new genomes have been copied to `Efs_assemblies/` (total: 8,603 genomes).

## Current Status (2025-12-19)

### All Pipeline Steps Complete ✅

- [x] New genomes copied to `Efs_assemblies/` (8,603 total)
- [x] Metadata file paths updated to `202251215_metadata_8573strains_23Isolates.txt`
- [x] **Step 01: Prokka annotation** - 8,603 genomes, avg 2,821 genes
- [x] **Step 03: BLAST search** - 96.2% with complete operon
- [x] **Step 04: Core gene analysis** - 1,251 core genes (≥95%)
- [x] **Step 05: Operon assembly extraction** - 7 genes, 8,275 avg seqs/gene
- [x] **Step 06: Comparative diversity analysis** - operon genes 6% more conserved
- [x] **Step 07: dN/dS analysis** - completed Dec 18
- [x] **Step 07b: Codon variation analysis** - completed Dec 19
- [x] **Step 07c: Visualize variation analysis** - completed Dec 18
- [x] **Step 08: Start site analysis** - completed Dec 18
- [x] **Step 09: Operon order analysis** - completed Dec 18
- [x] **Step 13e: Full operon mapping** - completed Dec 18
- [x] Git commits pushed for all updated files

### Remaining Tasks
- [ ] Draft manuscript results section using files in "Manuscript Statistics Files Reference" below
- [ ] Run verification script (`verify_analysis.sh`) for final check

### Issue Resolved: 21 Missing Genomes (Historical)

Job 9195589 initially processed only 8,582 of 8,603 genomes due to array task count bug.
**Fix applied**: Job 9198729 processed the 21 missing genomes.
**Final result**: All 8,603 genomes processed.

---

## Pipeline Steps (After Prokka Completes)

### Step 1: Verify Prokka Completion

```bash
# Check both jobs completed
sacct -j 9195589 --format=JobID,State,ExitCode | head -5
sacct -j 9198729 --format=JobID,State,ExitCode

# Count outputs (should be 8,603)
ls prokka_output | wc -l

# Optional: Validate outputs
cd 01_prokka_annotation
sbatch --array=1 run_prokka_pipeline.sh --validate
```

### Step 2: Regenerate Prokka Statistics

**IMPORTANT**: The initial stats file (`manuscript_stats_20251217_041626.json`) shows only 8,552 genomes because it was generated before the 21 missing genomes were processed.

**After job 9198729 completes, regenerate stats:**
```bash
cd 01_prokka_annotation
python prokka_pipeline.py stats
```

This creates a new `manuscript_stats_YYYYMMDD_HHMMSS.json` with all 8,603 genomes.

**Verify the new stats file:**
```bash
cd 01_prokka_annotation

# Find the latest stats file
ls -lt manuscript_stats_*.json | head -1

# View the stats
cat manuscript_stats_*.json | python -m json.tool
```

**Key statistics to extract for manuscript:**
- `input_count`: Total input assemblies
- `complete_genomes`: Successfully annotated genomes
- `avg_genes_per_genome`: Average genes per genome
- `total_genes`: Total genes predicted

**Note**: The old `manuscript_stats.txt` (from August 2025) is outdated. The `manuscript_numbers.py` script mentioned in the README is missing. You may need to manually format the JSON stats or recreate the script.

**Alternative - regenerate stats manually:**
```bash
cd 01_prokka_annotation
sbatch run_prokka_pipeline.sh --stats
```

### Step 3: BLAST Search (Step 03)

**Important**: Use `--prokka-dir ../prokka_output` to point to the new Prokka output location.

```bash
cd 03_blast_search

# Run BLAST search (array job, processes 100 genomes per batch)
sbatch --array=1-87 run_pipeline.sh --prokka-dir ../prokka_output
```

Note: Changed to `--array=1-87` because 8,601 genomes / 100 = 87 batches (was 86 for 8,587).

**BLAST modes that will run:**
- `coding_protein`: tblastn protein query vs genome assemblies
- `coding_nt`: blastn nucleotide query vs genome assemblies
- `prokka_variants`: blastn Prokka genes vs reference (uses `--prokka-dir`)
- `noncoding`: blastn promoter/regulatory vs genome assemblies

The pipeline **skips already processed genomes**, so only new ones will be BLASTed.

### Step 4: Process BLAST Results

After all BLAST array jobs complete:

```bash
cd 03_blast_search
sbatch run_pipeline.sh all
```

This runs:
- `process`: Combine BLAST results, identify best hits
- `overview`: Create gene prevalence statistics
- `stats`: Generate manuscript statistics

### Step 4: Core Gene Analysis (Step 04)

**Current status**: Results based on 8,587 genomes (old). Needs re-run for 8,601 genomes.

The core gene analysis identifies genes present in ≥95% of genomes as a baseline for comparing operon gene diversity.

**Input**: Prokka GFF files from `prokka_output/`
**Output**:
- `output/core_genes_95pct.txt` - List of core genes
- `output/gene_prevalence_stats.csv` - Prevalence stats for all genes
- `output/core_gene_sequences/` - FASTA files for each core gene

```bash
cd 04_core_gene_analysis

# Run with new Prokka output directory
sbatch run_pipeline.sh --prokka-dir ../prokka_output
```

**Expected changes**:
- Gene counts will change from 8,587 to 8,601
- Prevalence percentages will shift slightly (e.g., a gene in 8,158 genomes was 95.0%, now 94.8%)
- Some genes may cross the 95% threshold in either direction

**After completion, generate manuscript stats:**
```bash
python manuscript_numbers.py > manuscript_stats.txt
```

### Step 5: Operon Assembly Extraction (Step 05)

**Current status**: Results from August 22, 2025 with ~8,260 sequences per gene. Needs re-run for new genomes.

This step extracts operon gene sequences from genome assemblies based on BLAST coordinates, creates MSAs, and calculates conservation metrics.

**Input**:
- BLAST results: `../03_blast_search/output/blast_results/` (must be updated first!)
- Assemblies: `../../Efs_assemblies/` (already has 8,601 genomes)

**Current output (will change)**:
```
frpC: 8261 sequences → ~8,280+ expected
glpC: 8260 sequences → ~8,279+ expected
ptsD: 8261 sequences → ~8,280+ expected
ptsC: 8257 sequences → ~8,276+ expected
ptsB: 8260 sequences → ~8,279+ expected
ptsA: 8260 sequences → ~8,279+ expected
fruR: 8254 sequences → ~8,273+ expected
Promoter: 8586 sequences → ~8,600+ expected
```

```bash
cd 05_operon_assembly_extraction

# Run full pipeline (after BLAST completes)
sbatch run_operon_extraction.sh

# Or with gene boundary validation (uses Prokka annotations)
sbatch run_operon_extraction.sh --with-gene-boundary
```

**Note**: If using `--with-gene-boundary`, the script uses the old Prokka path (`../01_prokka_annotation/output/prokka_results`). The old Prokka outputs should still work for validation purposes, or you may need to update the script to use `../prokka_output`.

**After completion:**
```bash
# Check summary
cat output/extraction_pipeline_summary.txt

# Generate manuscript statistics
python manuscript_numbers.py > manuscript_stats.txt
```

### Step 6: Comparative Diversity Analysis (Step 06)

**Current status**: Uses pre-computed metrics from Steps 04 and 05. Will auto-update when dependencies re-run.

This step creates conservation ranking visualizations comparing operon genes to core genes.

**Input** (reads automatically):
- Core gene metrics: `../04_core_gene_analysis/output/core_gene_conservation_metrics.csv`
- Operon metrics: `../05_operon_assembly_extraction/output/operon_conservation_metrics_detailed.csv`

**Output**:
- `output/conservation_ranking_bars.png/pdf` - Bar plot showing top 20, operon genes, bottom 20
- `output/conservation_distribution.png/pdf` - Histogram with operon genes highlighted
- `output/operon_conservation_summary.csv` - Rankings table

**Current values** (will change after re-run):
- Total core genes: 1,252
- ptsA: Rank 153 (87.9th percentile) - most conserved operon gene
- glpC: Rank 815 (35.0th percentile) - least conserved operon gene

```bash
cd 06_diversity_analysis

# Re-run after Steps 04 and 05 complete
sbatch run_analysis.sh

# Or run directly
python comparative_analysis.py
```

**Files to check after re-run:**
```bash
cat output/operon_conservation_summary.csv  # Should show updated Total Core Genes count
```

### Step 7: dN/dS Analysis (Step 07)

**Current status**: Results in `output_old/` from August 2025. Needs re-run for 8,601 genomes.

This step calculates dN/dS ratios and selection metrics for operon and core genes.

**Input**:
- Operon alignments from Step 05: `../05_operon_assembly_extraction/output/msa/dna_alignments`
- Core gene alignments from Step 04: `../04_core_gene_analysis/output/core_gene_alignments`

**Output**:
- `output/tables/operon_dnds_analysis.csv` - dN/dS metrics for operon genes
- `output/tables/core_genes_dnds_analysis.csv` - dN/dS metrics for core genes
- `output/plots/` - Various dN/dS plots

```bash
cd 07_dnds_analysis

# Run after Steps 04 and 05 complete
sbatch run_dnds_analysis.sh
```

**Note**: This is a long-running job (72h time limit, 64GB memory).

### Step 7b: Codon Variation Analysis (Step 07b)

**Current status**: Generates `pN_pS_scatter.pdf`. Depends on Step 07 output.

**Input** (from Step 07):
- `../07_dnds_analysis/output/tables/operon_dnds_analysis.csv`
- `../07_dnds_analysis/output/tables/core_genes_dnds_analysis.csv`

**Output**:
- `output/plots/pN_pS_scatter.pdf` - Scatter plot of pN vs pS
- `output/plots/selection_distribution.pdf` - Distribution comparison
- `output/codon_variation_report.txt` - Summary report

```bash
cd 07b_codon_variation_analysis

# Run after Step 07 completes
sbatch run_codon_variation.sh
```

**Path fix applied**: Updated input paths from `output_old/` to `output/` so it picks up new Step 07 results.

### Step 13e: Full Operon Mapping

**Current status**: Results from October 2025 with ~8,587 genomes. Needs re-run for 8,601 genomes.

This step maps the full operon sequence to each assembly and compares legacy vs SNP-aware reference.

**Input**:
- Assemblies: `../../Efs_assemblies` (will automatically scan all 8,601 genomes)
- Legacy GenBank: `../02_reference_operon_extraction/operon.gb` (reference, no change)
- Updated GenBank: `../13a_new_reference_operon_extraction/FL_operon_with_SNPs.gb` (reference, no change)

**Run order** (must be sequential):
```bash
cd 13e_full_operon_mapping

# 1. Main mapping (creates output_full_run/)
sbatch run_full_mapping.sh

# 2. Process raw BLAST output (after step 1 completes)
sbatch run_process_raw.sh

# 3. Generate position plots (after step 2 completes)
sbatch run_position_plots.sh

# 4. Analyze variants (after step 1 completes)
sbatch run_analyze_variants.sh

# 5. Generate manuscript statistics (after all above complete)
python manuscript_numbers.py
```

**Output files to check**:
```
output_full_run/
├── full_operon_mapping.tsv          # Per-assembly hit metrics
├── summary.txt                       # Aggregated hit counts
├── legacy_position_summary.tsv       # Legacy position-level stats
├── updated_position_summary.tsv      # Updated position-level stats
├── variant_site_summary.tsv          # Variant site statistics
├── variant_site_allele_counts.tsv    # Per-site allele counts
├── variant_site_calls.tsv            # Per-assembly allele calls
├── manuscript_numbers.txt            # Final manuscript statistics
├── position_stats/
│   ├── position_metric_summary.tsv   # Combined metrics
│   └── plots/
│       ├── legacy_coverage.png
│       ├── legacy_mismatch_rate.png
│       ├── legacy_deletion_rate.png
│       ├── updated_coverage.png
│       ├── updated_mismatch_rate.png
│       └── updated_deletion_rate.png
└── raw_blast/
    ├── legacy/*.tsv                  # Raw BLAST output per assembly
    └── updated/*.tsv
```

**Key statistics to verify after re-run**:
- `summary.txt`: Should show ~8,601 assemblies scanned
- `manuscript_numbers.txt`: Should report updated counts

### Step 7+: Additional Downstream Steps

After Steps 03-06 complete, these steps also need re-running:

| Step | Directory | Depends On | Re-run Command | Status |
|------|-----------|------------|----------------|--------|
| 07 | `07_dnds_analysis/` | Steps 04 + 05 | `sbatch run_dnds_analysis.sh` | ✅ Paths OK |
| 07b | `07b_codon_variation_analysis/` | Step 07 | `sbatch run_codon_variation.sh` | ✅ Fixed (was output_old/) |
| 07c | `07c_visualize_variation_analysis/` | Step 05 MSAs | `sbatch run_visualization.sh --genes ptsA ptsB ptsC ptsD frpC glpC fruR` | ✅ Paths OK |
| 08 | `08_start_site_analysis/` | Prokka + metadata | `sbatch run_start_site_analysis.sh` | ✅ Fixed (was old Prokka path) |
| 09 | `09_operon_order_analysis/` | BLAST results | `sbatch run_pipeline.sh` | ✅ Paths OK |
| 13e | `13e_full_operon_mapping/` | BLAST + assemblies | `sbatch run_full_mapping.sh` | ✅ Paths OK |

**Path verification summary:**
- **07**: Uses `../05_operon_assembly_extraction/output/msa/dna_alignments` and `../04_core_gene_analysis/output/core_gene_alignments` ✅
- **07b**: Updated from `output_old/` to `output/` ✅
- **07c**: Uses `../05_operon_assembly_extraction/output/msa/dna_alignments` ✅
- **08**: Updated `PROKKA_DIR` from old path to `../prokka_output` ✅
- **09**: Uses `../03_blast_search/output/blast_results` ✅
- **13e**: Uses `../../Efs_assemblies` (assemblies directly) ✅

**Dependency chain:**
```
Prokka (01) ─┬─> BLAST (03) ─┬─> Operon Extraction (05) ─┬─> Comparative (06)
             │               │                           ├─> Variation (07c)
             │               │                           └─> dN/dS (07) ──> Codon Var (07b)
             │               └─> Operon Order (09)
             │               └─> Full Mapping (13e)
             │
             ├─> Core Genes (04) ──┬──────────────────────> Comparative (06)
             │                     └──────────────────────> dN/dS (07)
             │
             └─> Start Site (08)
```

---

## Key Findings

### Path Discrepancy
There are **two** Prokka output locations:
- `01_prokka_annotation/output/prokka_results/` - Old results (July 2025), 8,589 genomes
- `prokka_output/` - New results (current job), will have 8,601 genomes

The BLAST pipeline defaults to the old location. **Always use `--prokka-dir ../prokka_output`** for the new analysis.

### Old Manuscript Stats
The file `01_prokka_annotation/manuscript_stats.txt` is from August 2025 and shows 8,587 genomes. This will be regenerated after Prokka completes.

### Metadata File
Updated from `8587_Efs_metadata_ASbarcode.txt` to `202251215_metadata_8573strains_23Isolates.txt` in all scripts.

---

## Quick Reference Commands

```bash
# Check Prokka job status
squeue -j 9195589

# Count completed Prokka outputs
ls prokka_output | wc -l

# After Prokka done - run BLAST
cd 03_blast_search && sbatch --array=1-87 run_pipeline.sh --prokka-dir ../prokka_output

# After BLAST done - process results
cd 03_blast_search && sbatch run_pipeline.sh all

# After BLAST done - run Core Gene Analysis
cd 04_core_gene_analysis && sbatch run_pipeline.sh --prokka-dir ../prokka_output

# After BLAST done - run Operon Extraction
cd 05_operon_assembly_extraction && sbatch run_operon_extraction.sh

# After Steps 04 and 05 - run Comparative Analysis
cd 06_diversity_analysis && sbatch run_analysis.sh

# After Steps 04 and 05 - run dN/dS Analysis (long-running job)
cd 07_dnds_analysis && sbatch run_dnds_analysis.sh

# After Step 07 - run Codon Variation Analysis (generates pN_pS_scatter.pdf)
cd 07b_codon_variation_analysis && sbatch run_codon_variation.sh

# After Step 05 - run downstream steps
cd 07c_visualize_variation_analysis && sbatch run_visualization.sh --genes ptsA ptsB ptsC ptsD frpC glpC fruR
cd 08_start_site_analysis && sbatch run_start_site_analysis.sh
cd 09_operon_order_analysis && sbatch run_pipeline.sh

# Step 13e - Full Operon Mapping (run in order, wait for each to complete)
cd 13e_full_operon_mapping
sbatch run_full_mapping.sh          # Wait for completion, then:
sbatch run_process_raw.sh           # Wait for completion, then:
sbatch run_position_plots.sh        # Can run in parallel with next:
sbatch run_analyze_variants.sh      # After all complete:
python manuscript_numbers.py        # Generate final stats

# Check any job status
squeue -u $USER

# After ALL steps complete - run verification
cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis
./verify_analysis.sh 2025-12-16
```

## Statistics to Collect for Manuscript

After all steps complete, collect these statistics:

### Files to Check After Each Step

| Step | Output Files | Regenerate Stats Command |
|------|--------------|--------------------------|
| 01 Prokka | `01_prokka_annotation/manuscript_stats_*.json` | `sbatch run_prokka_pipeline.sh --stats` |
| 03 BLAST | `03_blast_search/output/manuscript_stats.txt`<br>`03_blast_search/output/operon_summary.csv`<br>`03_blast_search/output/gene_prevalence.csv` | `sbatch run_pipeline.sh stats` |
| 04 Core Genes | `04_core_gene_analysis/output/gene_prevalence_stats.csv`<br>`04_core_gene_analysis/output/core_genes_95pct.txt` | `python manuscript_numbers.py > manuscript_stats.txt` |
| 05 Operon Extraction | `05_operon_assembly_extraction/output/extraction_pipeline_summary.txt`<br>`05_operon_assembly_extraction/output/operon_conservation_metrics.csv`<br>`05_operon_assembly_extraction/manuscript_stats.txt` | `python manuscript_numbers.py > manuscript_stats.txt` |
| 06 Comparative Analysis | `06_diversity_analysis/output/operon_conservation_summary.csv`<br>`06_diversity_analysis/output/conservation_ranking_bars.pdf` | `python comparative_analysis.py` |
| 07 dN/dS | `07_dnds_analysis/output/tables/operon_dnds_analysis.csv`<br>`07_dnds_analysis/output/tables/core_genes_dnds_analysis.csv` | `sbatch run_dnds_analysis.sh` |
| 07b Codon Variation | `07b_codon_variation_analysis/output/plots/pN_pS_scatter.pdf`<br>`07b_codon_variation_analysis/output/codon_variation_report.txt` | `sbatch run_codon_variation.sh` |
| 08 Start Site | `08_start_site_analysis/output/start_site_summary.tsv`<br>`08_start_site_analysis/output/manuscript_stats.txt`<br>`08_start_site_analysis/output/metadata_stratification_summary.txt` | `python manuscript_numbers.py` |
| 13e Full Mapping | `13e_full_operon_mapping/output_full_run/summary.txt`<br>`13e_full_operon_mapping/output_full_run/full_operon_mapping.tsv`<br>`13e_full_operon_mapping/output_full_run/manuscript_numbers.txt`<br>`13e_full_operon_mapping/output_full_run/position_stats/plots/*.png` | `python manuscript_numbers.py` |

### Key Numbers to Update in Manuscript

**From Step 01 (Prokka):**
- Total genomes: 8,587 → **8,601** (change: +14)
- Average genes per genome
- Total genes predicted

**From Step 03 (BLAST):**
- Genomes with complete operon (7 genes)
- Gene prevalence percentages

**From Step 04 (Core Genes):**
- Number of core genes at 95% threshold
- Gene counts will change from 8,587 → 8,601

**From Step 05 (Operon Extraction):**
Actual values (Dec 2025):
- Sequences per gene: ~8,275 average
- ptsA genes analyzed: 8,167
- Laboratory strains: n=324
- Other niches: n=5,703 (with niche data)

**Key Finding VERIFIED (Dec 2025):**
The ptsA laboratory adaptation pattern (TTG→ATG start codon) remains significant:
- Old (8,587 genomes): Laboratory ATG 55.1% vs Other 10.9% (p=4.4e-79)
- New (8,603 genomes): Laboratory ATG 60.5% vs Other 21.3% (p=3.15e-50, OR=5.86)
- ✅ Pattern confirmed: Lab strains 5.9x more likely to use ATG start codon

### Final Verification Checklist

After ALL steps complete, run this verification script to ensure everything worked:

```bash
#!/bin/bash
# Save as verify_analysis.sh and run from operon_analysis/ directory
# Usage: ./verify_analysis.sh [YYYY-MM-DD]
# If date provided, checks that files were modified AFTER that date

echo "=============================================="
echo "ANALYSIS VERIFICATION CHECKLIST"
echo "=============================================="

# Reference date for checking if files are updated
# Default: 2025-12-16 (when re-run started)
REF_DATE="${1:-2025-12-16}"
echo "Checking files modified after: $REF_DATE"
echo ""

# Function to check file exists and is newer than reference date
check_file() {
    local file="$1"
    local desc="$2"

    if [ ! -f "$file" ]; then
        echo "  ❌ MISSING: $desc"
        return 1
    fi

    # Get file modification date
    FILE_DATE=$(stat -c %Y "$file" 2>/dev/null)
    REF_EPOCH=$(date -d "$REF_DATE" +%s 2>/dev/null)

    if [ -n "$FILE_DATE" ] && [ -n "$REF_EPOCH" ] && [ "$FILE_DATE" -lt "$REF_EPOCH" ]; then
        MOD_DATE=$(stat -c %y "$file" | cut -d' ' -f1)
        echo "  ⚠️  STALE: $desc (last modified: $MOD_DATE)"
        return 2
    else
        MOD_DATE=$(stat -c %y "$file" | cut -d' ' -f1)
        echo "  ✅ OK: $desc ($MOD_DATE)"
        return 0
    fi
}

# Step 01: Prokka
echo "=== Step 01: Prokka Annotation ==="
PROKKA_COUNT=$(ls prokka_output 2>/dev/null | wc -l)
echo "  Prokka outputs: $PROKKA_COUNT (expected: ~8,601)"
LATEST_JSON=$(ls -t 01_prokka_annotation/manuscript_stats_*.json 2>/dev/null | head -1)
if [ -n "$LATEST_JSON" ]; then
    check_file "$LATEST_JSON" "manuscript_stats JSON"
else
    echo "  ❌ MISSING: No manuscript_stats JSON found"
fi
echo ""

# Step 03: BLAST
echo "=== Step 03: BLAST Search ==="
check_file "03_blast_search/output/manuscript_stats.txt" "manuscript_stats.txt"
check_file "03_blast_search/output/operon_summary.csv" "operon_summary.csv"
check_file "03_blast_search/output/gene_prevalence.csv" "gene_prevalence.csv"
echo ""

# Step 04: Core Genes
echo "=== Step 04: Core Gene Analysis ==="
check_file "04_core_gene_analysis/output/gene_prevalence_stats.csv" "gene_prevalence_stats.csv"
check_file "04_core_gene_analysis/output/core_genes_95pct.txt" "core_genes_95pct.txt"
check_file "04_core_gene_analysis/output/core_gene_conservation_metrics.csv" "core_gene_conservation_metrics.csv"
if [ -f "04_core_gene_analysis/output/gene_prevalence_stats.csv" ]; then
    CORE_COUNT=$(head -2 04_core_gene_analysis/output/gene_prevalence_stats.csv | tail -1 | cut -d',' -f3)
    echo "  📊 Genome count in stats: $CORE_COUNT (expected: 8,601)"
fi
echo ""

# Step 05: Operon Extraction
echo "=== Step 05: Operon Extraction ==="
check_file "05_operon_assembly_extraction/output/extraction_pipeline_summary.txt" "extraction_pipeline_summary.txt"
check_file "05_operon_assembly_extraction/output/operon_conservation_metrics.csv" "operon_conservation_metrics.csv"
check_file "05_operon_assembly_extraction/manuscript_stats.txt" "manuscript_stats.txt"
if [ -f "05_operon_assembly_extraction/output/extraction_pipeline_summary.txt" ]; then
    echo "  📊 Sequences extracted:"
    grep -E "^\s+(frpC|ptsA|ptsD):" 05_operon_assembly_extraction/output/extraction_pipeline_summary.txt 2>/dev/null | head -3 | sed 's/^/      /'
fi
echo ""

# Step 06: Comparative Analysis
echo "=== Step 06: Comparative Analysis ==="
check_file "06_diversity_analysis/output/operon_conservation_summary.csv" "operon_conservation_summary.csv"
check_file "06_diversity_analysis/output/conservation_ranking_bars.pdf" "conservation_ranking_bars.pdf"
if [ -f "06_diversity_analysis/output/operon_conservation_summary.csv" ]; then
    TOTAL_CORE=$(head -2 06_diversity_analysis/output/operon_conservation_summary.csv | tail -1 | cut -d',' -f4)
    echo "  📊 Total core genes: $TOTAL_CORE"
fi
echo ""

# Step 07: dN/dS Analysis
echo "=== Step 07: dN/dS Analysis ==="
check_file "07_dnds_analysis/output/tables/operon_dnds_analysis.csv" "operon_dnds_analysis.csv"
check_file "07_dnds_analysis/output/tables/core_genes_dnds_analysis.csv" "core_genes_dnds_analysis.csv"
echo ""

# Step 07b: Codon Variation Analysis
echo "=== Step 07b: Codon Variation Analysis ==="
check_file "07b_codon_variation_analysis/output/plots/pN_pS_scatter.pdf" "pN_pS_scatter.pdf"
check_file "07b_codon_variation_analysis/output/plots/selection_distribution.pdf" "selection_distribution.pdf"
check_file "07b_codon_variation_analysis/output/codon_variation_report.txt" "codon_variation_report.txt"
echo ""

# Step 07c: Variation Visualization
echo "=== Step 07c: Variation Visualization ==="
PLOT_COUNT=$(ls 07c_visualize_variation_analysis/output/plots/*_variable_codons.pdf 2>/dev/null | wc -l)
TABLE_COUNT=$(ls 07c_visualize_variation_analysis/output/tables/*_variable_codons.csv 2>/dev/null | wc -l)
echo "  📊 Variable codon plots: $PLOT_COUNT files"
echo "  📊 Variable codon tables: $TABLE_COUNT files"
check_file "07c_visualize_variation_analysis/output/plots/ptsA_variable_codons.pdf" "ptsA_variable_codons.pdf"
echo ""

# Step 08: Start Site Analysis
echo "=== Step 08: Start Site Analysis ==="
check_file "08_start_site_analysis/output/start_site_summary.tsv" "start_site_summary.tsv"
check_file "08_start_site_analysis/output/manuscript_stats.txt" "manuscript_stats.txt"
check_file "08_start_site_analysis/output/metadata_stratification_summary.txt" "metadata_stratification_summary.txt"
if [ -f "08_start_site_analysis/output/start_site_summary.tsv" ]; then
    PTSA_COUNT=$(grep -c "ptsA" 08_start_site_analysis/output/start_site_summary.tsv 2>/dev/null || echo "0")
    echo "  📊 ptsA entries: $PTSA_COUNT (expected: ~8,280+)"
fi
echo ""

# Step 09: Operon Order
echo "=== Step 09: Operon Order Analysis ==="
check_file "09_operon_order_analysis/output/operon_order_summary.csv" "operon_order_summary.csv"
check_file "09_operon_order_analysis/output/manuscript_stats.txt" "manuscript_stats.txt"
echo ""

# Step 13e: Full Operon Mapping
echo "=== Step 13e: Full Operon Mapping ==="
MAPPING_CSV=$(ls 13e_full_operon_mapping/output*/mapping_summary.csv 2>/dev/null | head -1)
MAPPING_TXT=$(ls 13e_full_operon_mapping/output*/manuscript_numbers.txt 2>/dev/null | head -1)
[ -n "$MAPPING_CSV" ] && check_file "$MAPPING_CSV" "mapping_summary.csv" || echo "  ❌ MISSING: mapping_summary.csv"
[ -n "$MAPPING_TXT" ] && check_file "$MAPPING_TXT" "manuscript_numbers.txt" || echo "  ❌ MISSING: manuscript_numbers.txt"
echo ""

echo "=============================================="
echo "VERIFICATION COMPLETE"
echo "=============================================="
echo ""
echo "Legend: ✅ OK (updated) | ⚠️ STALE (old file) | ❌ MISSING"
echo ""
echo "Next steps:"
echo "  - Re-run any steps with STALE or MISSING files"
echo "  - Regenerate manuscript_stats after fixing issues"
```

### Quick Verification Commands

```bash
# Check timestamps to ensure files are fresh (should be after re-run date)
ls -la 01_prokka_annotation/manuscript_stats_*.json
ls -la 03_blast_search/output/manuscript_stats.txt
ls -la 04_core_gene_analysis/output/gene_prevalence_stats.csv
ls -la 05_operon_assembly_extraction/output/extraction_pipeline_summary.txt
ls -la 05_operon_assembly_extraction/manuscript_stats.txt
ls -la 06_diversity_analysis/output/operon_conservation_summary.csv
ls -la 07c_visualize_variation_analysis/output/plots/ptsA_variable_codons.pdf
ls -la 08_start_site_analysis/output/manuscript_stats.txt
ls -la 09_operon_order_analysis/output/manuscript_stats.txt
ls -la 13e_full_operon_mapping/output*/manuscript_numbers.txt

# Verify genome counts in key files
head -5 04_core_gene_analysis/output/gene_prevalence_stats.csv  # Should show count=8601
head -20 05_operon_assembly_extraction/output/extraction_pipeline_summary.txt
```

---

## Files Modified

| File | Change |
|------|--------|
| `JOBS.md` | Created - tracks running jobs |
| `PLAN.md` | Created - this file |
| `11_draw_tree/draw_tree_with_metadata.py` | Updated metadata path |
| `11_draw_tree/README.md` | Updated metadata path |
| `08_start_site_analysis/analyze_start_sites.py` | Updated metadata path (3 places) |
| `08_start_site_analysis/run_start_site_analysis.sh` | Updated metadata path |
| `08_start_site_analysis/manuscript_numbers.py` | Updated metadata path |
| `09_operon_order_analysis/operon_order_pipeline.py` | Updated metadata path |
| `12_visualize_mapping/generate_visualization.py` | Updated metadata path |
| `05_operon_assembly_extraction/manuscript_numbers.py` | Updated metadata path |
| `13c_new_operon_assembly_extraction/manuscript_numbers.py` | Updated metadata path |
| `.gitignore` | Updated metadata path |
| `08_start_site_analysis/run_start_site_analysis.sh` | Fixed PROKKA_DIR path (`../prokka_output`) |
| `07b_codon_variation_analysis/run_codon_variation.sh` | Fixed input paths from `output_old/` to `output/` |
| `verify_analysis.sh` | Created - verification script to check all outputs |
| `verify_analysis.sh` | Updated - added Steps 07 and 07b checks |

---

## Manuscript Statistics Files Reference (Dec 2025)

Use these files to draft the analysis results section. All files updated Dec 17-19, 2025 with 8,603 genomes.

### Primary Statistics Files

| Step | File | Updated | Description |
|------|------|---------|-------------|
| **01 Prokka** | `01_prokka_annotation/manuscript_stats.txt` | Dec 17 | 8,603 genomes, avg 2,821 genes, median 2,800, range 2,391-3,728, SD 201 |
| **03 BLAST** | `03_blast_search/output/manuscript_stats.txt` | Dec 17 | Operon gene prevalence across genomes |
| **04 Core Genes** | `04_core_gene_analysis/manuscript_stats.txt` | Dec 18 | 1,251 core genes at 95% threshold |
| **05 Operon Extraction** | `05_operon_assembly_extraction/manuscript_stats.txt` | Dec 18 | 7 genes, ~8,275 sequences/gene |
| **06 Diversity** | `06_diversity_analysis/manuscript_stats.txt` | Dec 18 | Conservation rankings, operon 6% more conserved |
| **07b Codon Variation** | `07b_codon_variation_analysis/output/codon_variation_report.txt` | Dec 19 | pN/pS selection analysis |
| **08 Start Sites** | `08_start_site_analysis/output/manuscript_stats.txt` | Dec 18 | Start codon usage by niche |
| **09 Operon Order** | `09_operon_order_analysis/output/manuscript_stats.txt` | Dec 18 | Gene order patterns |
| **13e Full Mapping** | `13e_full_operon_mapping/output_full_run/manuscript_numbers.txt` | Dec 18 | Position-level mapping stats |

### Supporting Data Files

| File | Updated | Content |
|------|---------|---------|
| `03_blast_search/output/gene_prevalence.csv` | Dec 17 | Per-gene presence counts |
| `03_blast_search/output/operon_summary.csv` | Dec 17 | Per-genome operon completeness |
| `04_core_gene_analysis/output/core_gene_conservation_summary.txt` | Dec 17 | Core gene metrics |
| `05_operon_assembly_extraction/output/operon_conservation_metrics_detailed.csv` | Dec 18 | Per-gene conservation |
| `06_diversity_analysis/output/operon_conservation_summary.csv` | Dec 18 | Rankings table |
| `08_start_site_analysis/output/metadata_stratification_summary.txt` | Dec 18 | Sample counts by niche |

### Python Scripts (to regenerate stats)

| Step | Script |
|------|--------|
| 04 | `04_core_gene_analysis/manuscript_numbers.py` |
| 05 | `05_operon_assembly_extraction/manuscript_numbers.py` |
| 06 | `06_diversity_analysis/manuscript_numbers.py` |
| 08 | `08_start_site_analysis/manuscript_numbers.py` |
| 13e | `13e_full_operon_mapping/manuscript_numbers.py` |

### Old Files (Ignore)

- `02_reference_operon_extraction/manuscript_stats.txt` (Aug 21) - reference operon only
- `05_operon_assembly_extraction/manuscript_numbers.txt` (Aug 22) - superseded by .py
- `01_prokka_annotation/manuscript_stats_2025111*.json` (Nov) - old test runs
