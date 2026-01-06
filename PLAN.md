# Analysis Re-run Plan: Correct Genome Count (8,596)

## Current Progress (2026-01-06)

| Step | Status | Result |
|------|--------|--------|
| **1.1 Prokka Stats** | ✅ Complete | 8,596 genomes, avg 2820.8 genes |
| **1.2 BLAST Results** | ✅ Complete | operon_summary.csv filtered (8,596 genomes) |
| **1.2b BLAST CSVs** | ⏳ Pending | all_blast_hits*.csv, pts*_filtered_hits.csv need filtering |
| **1.3 Core Gene MSAs** | ✅ Complete | 1,251 files filtered, 8,020 seqs removed, conservation=0.9065 |
| **1.4 Operon MSAs** | ✅ Complete | 7 files filtered, 49 seqs removed, conservation=0.9878 |
| **2.1 Diversity Analysis** | 🔄 Running | SLURM job 9232732 |
| **2.2 dN/dS Analysis** | ⏳ Pending | (slow - consider skipping) |
| **2.5 Start Site Analysis** | ⏳ Pending | |
| **2.6 Operon Order Analysis** | ⏳ Pending | |

**Phase 1 mostly complete.** Diversity analysis running. BLAST detail CSVs still need filtering.

### Pending tasks:
- `sbatch run_filter_remaining.sh` - filter BLAST detail CSVs
- `git push origin main` - push 3 commits to GitHub

---

## Overview

**ISSUE IDENTIFIED (2026-01-06)**: The analysis includes 7 extra assemblies not in metadata.

| Current | Correct | Difference |
|---------|---------|------------|
| 8,603   | 8,596   | -7 (0.08%) |

## Key Insight: Efficient Re-run Strategy

**We do NOT need to re-run everything from scratch!**

- **Prokka**: Already done - just filter stats
- **BLAST**: Already done - just filter output CSVs
- **MSAs**: Can be FILTERED (remove 7 sequences) without re-alignment
- **Downstream**: Re-run with filtered inputs

This saves hours of compute time (especially avoiding re-running MAFFT on 1,251+ alignments).

---

## Assemblies to EXCLUDE (7 total)

```
ENT_CA1913AA_AS.genome
ENT_CA1914AA_AS.genome
ENT_CA1915AA_AS.genome
ENT_CA1916AA_AS.genome
ENT_CA1917AA_AS.genome
ENT_CA1918AA_AS.genome
ENT_CA1919AA_AS.genome
```

Helper files created:
- `exclude_assemblies.txt` - List of 7 to exclude
- `valid_as_barcodes.txt` - List of 8,596 valid IDs

---

## Step-by-Step Re-run Plan

### Phase 1: Filter Existing Outputs (No Re-computation)

#### Step 1.1: Prokka Stats (01_prokka_annotation)

**Action**: Regenerate stats by counting only valid genomes

```bash
cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/01_prokka_annotation
conda activate efs_diversity

python3 << 'PYTHON'
import os
import statistics
from datetime import datetime

# Load valid AS_Barcodes
with open('../valid_as_barcodes.txt') as f:
    valid_ids = set(line.strip() for line in f)

print(f"Valid AS_Barcodes: {len(valid_ids)}")

# Count genes from valid prokka outputs only
prokka_dir = '../prokka_output'
gene_counts = []
skipped = []

for dirname in os.listdir(prokka_dir):
    # Normalize name (remove .result or .genome suffix)
    base_id = dirname.replace('.result', '').replace('.genome', '')

    if base_id not in valid_ids:
        skipped.append(dirname)
        continue

    # Count CDS from GFF file
    gff_path = os.path.join(prokka_dir, dirname, f"{dirname}.gff")
    if os.path.exists(gff_path):
        gene_count = sum(1 for line in open(gff_path) if '\tCDS\t' in line)
        gene_counts.append(gene_count)

print(f"Processed: {len(gene_counts)} genomes")
print(f"Skipped: {len(skipped)} ({skipped})")

stats_dict = {
    'genomes': len(gene_counts),
    'avg': round(statistics.mean(gene_counts), 1),
    'median': int(statistics.median(gene_counts)),
    'min': min(gene_counts),
    'max': max(gene_counts),
    'stdev': round(statistics.stdev(gene_counts), 1),
    'total': sum(gene_counts)
}

# Write updated stats
with open('manuscript_stats.txt', 'w') as f:
    f.write(f"Prokka Annotation Statistics\n")
    f.write(f"============================\n")
    f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    f.write(f"Genomes analyzed: {stats_dict['genomes']:,}\n")
    f.write(f"Average genes per genome: {stats_dict['avg']}\n")
    f.write(f"Median genes per genome: {stats_dict['median']}\n")
    f.write(f"Gene count range: {stats_dict['min']} - {stats_dict['max']}\n")
    f.write(f"Standard deviation: {stats_dict['stdev']}\n")
    f.write(f"Total genes predicted: {stats_dict['total']:,}\n")

print(f"\nWrote manuscript_stats.txt")
print(f"Genomes: {stats_dict['genomes']:,}")
print(f"Avg genes: {stats_dict['avg']}")
PYTHON
```

---

#### Step 1.2: BLAST Results (03_blast_search)

**Action**: Filter operon_summary.csv and gene_prevalence.csv

```bash
cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/03_blast_search
conda activate efs_diversity

python3 << 'PYTHON'
import pandas as pd
from datetime import datetime

# Load valid AS_Barcodes
with open('../valid_as_barcodes.txt') as f:
    valid_ids = set(line.strip() for line in f)

print(f"Filtering to {len(valid_ids)} valid genomes")

# Read operon_summary.csv
df = pd.read_csv('output/operon_summary.csv')
print(f"Original rows: {len(df)}")

# Normalize genome IDs and filter
df['genome_id'] = df['genome'].str.replace('.result', '', regex=False).str.replace('.genome', '', regex=False)
df_filtered = df[df['genome_id'].isin(valid_ids)].copy()
print(f"Filtered rows: {len(df_filtered)}")
print(f"Removed: {len(df) - len(df_filtered)}")

# Save filtered (drop helper column)
df_filtered.drop(columns=['genome_id']).to_csv('output/operon_summary.csv', index=False)

# Recalculate gene_prevalence.csv
gene_cols = [c for c in df_filtered.columns if c not in ['genome', 'genome_id', 'total_genes', 'has_complete_operon']]
prevalence = []
for gene in gene_cols:
    if df_filtered[gene].dtype in ['int64', 'float64', 'bool']:
        count = int(df_filtered[gene].sum())
        pct = round(100 * count / len(df_filtered), 2)
        prevalence.append({'gene': gene, 'count': count, 'percentage': pct})

prev_df = pd.DataFrame(prevalence)
prev_df.to_csv('output/gene_prevalence.csv', index=False)
print(f"\nGene prevalence updated:")
print(prev_df.to_string(index=False))

# Update manuscript_stats.txt
total = len(df_filtered)
complete = int(df_filtered['has_complete_operon'].sum()) if 'has_complete_operon' in df_filtered.columns else 0

with open('output/manuscript_stats.txt', 'w') as f:
    f.write(f"BLAST Search Statistics\n")
    f.write(f"=======================\n")
    f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    f.write(f"Total genomes analyzed: {total:,}\n")
    f.write(f"Genomes with complete operon: {complete:,} ({100*complete/total:.1f}%)\n\n")
    f.write(f"Gene prevalence:\n")
    for _, row in prev_df.iterrows():
        f.write(f"  {row['gene']}: {row['count']:,} ({row['percentage']:.1f}%)\n")

print(f"\nUpdated manuscript_stats.txt")
PYTHON
```

---

#### Step 1.3: Filter Core Gene MSAs (04_core_gene_analysis)

**Action**: Remove excluded sequences from all 1,251 MSA files, then recalculate conservation

```bash
cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/04_core_gene_analysis
conda activate efs_diversity

python3 << 'PYTHON'
import os
from Bio import SeqIO
from datetime import datetime
import pandas as pd
import numpy as np
from collections import Counter

# Load exclusion list (normalized IDs)
exclude_patterns = [
    'ENT_CA1913AA_AS', 'ENT_CA1914AA_AS', 'ENT_CA1915AA_AS',
    'ENT_CA1916AA_AS', 'ENT_CA1917AA_AS', 'ENT_CA1918AA_AS', 'ENT_CA1919AA_AS'
]

def should_exclude(seq_id):
    """Check if sequence ID matches any exclusion pattern"""
    for pattern in exclude_patterns:
        if pattern in seq_id:
            return True
    return False

def calculate_conservation(alignment_file):
    """Calculate Shannon entropy-based conservation score"""
    sequences = list(SeqIO.parse(alignment_file, 'fasta'))
    if not sequences:
        return None, 0, 0

    aln_length = len(sequences[0].seq)
    n_seqs = len(sequences)

    # Calculate per-position entropy
    entropies = []
    for pos in range(aln_length):
        column = [str(seq.seq[pos]).upper() for seq in sequences]
        counts = Counter(column)
        total = sum(counts.values())

        # Shannon entropy
        entropy = 0
        for count in counts.values():
            if count > 0:
                p = count / total
                entropy -= p * np.log2(p)

        # Normalize (max entropy for 4 bases = 2 bits)
        max_entropy = np.log2(min(4, total))
        if max_entropy > 0:
            norm_entropy = entropy / max_entropy
        else:
            norm_entropy = 0
        entropies.append(norm_entropy)

    # Conservation = 1 - mean(normalized_entropy)
    mean_entropy = np.mean(entropies)
    conservation = 1 - mean_entropy

    return conservation, n_seqs, aln_length

# Process all MSA files
msa_dir = 'output/core_gene_alignments'
msa_files = [f for f in os.listdir(msa_dir) if f.endswith('.fasta')]
print(f"Processing {len(msa_files)} MSA files...")

results = []
filtered_count = 0

for i, msa_file in enumerate(msa_files):
    if (i + 1) % 100 == 0:
        print(f"  {i + 1}/{len(msa_files)}...")

    filepath = os.path.join(msa_dir, msa_file)
    gene_name = msa_file.replace('.fasta', '')

    # Read sequences, filter out excluded ones
    sequences = list(SeqIO.parse(filepath, 'fasta'))
    original_count = len(sequences)

    filtered_seqs = [seq for seq in sequences if not should_exclude(seq.id)]
    removed = original_count - len(filtered_seqs)

    if removed > 0:
        filtered_count += 1
        # Write back filtered alignment
        SeqIO.write(filtered_seqs, filepath, 'fasta')

    # Calculate conservation on filtered alignment
    conservation, n_seqs, aln_len = calculate_conservation(filepath)

    results.append({
        'gene': gene_name,
        'n_sequences': n_seqs,
        'alignment_length': aln_len,
        'conservation_score': conservation,
        'sequences_removed': removed
    })

print(f"\nFiltered {filtered_count} MSA files")

# Save updated conservation metrics
df = pd.DataFrame(results)
df = df.sort_values('conservation_score', ascending=False)
df.to_csv('output/core_gene_conservation_metrics.csv', index=False)

# Update summary
total_seqs_removed = df['sequences_removed'].sum()
print(f"Total sequences removed: {total_seqs_removed}")
print(f"Mean conservation: {df['conservation_score'].mean():.4f}")
print(f"Median sequences per gene: {df['n_sequences'].median():.0f}")

# Update manuscript_stats.txt
with open('manuscript_stats.txt', 'w') as f:
    f.write(f"Core Gene Analysis Statistics\n")
    f.write(f"=============================\n")
    f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    f.write(f"Genomes analyzed: 8,596\n")
    f.write(f"Core genes (≥95% prevalence): {len(df):,}\n")
    f.write(f"Mean conservation score: {df['conservation_score'].mean():.3f}\n")
    f.write(f"Average sequences per gene: {df['n_sequences'].mean():.0f}\n")
    f.write(f"\nNote: Filtered {total_seqs_removed} sequences from excluded assemblies\n")

print("\nDone! Updated manuscript_stats.txt and conservation_metrics.csv")
PYTHON
```

---

#### Step 1.4: Filter Operon MSAs (05_operon_assembly_extraction)

**Action**: Same approach - filter MSAs and recalculate conservation

```bash
cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/05_operon_assembly_extraction
conda activate efs_diversity

python3 << 'PYTHON'
import os
from Bio import SeqIO
from datetime import datetime
import pandas as pd
import numpy as np
from collections import Counter

exclude_patterns = [
    'ENT_CA1913AA_AS', 'ENT_CA1914AA_AS', 'ENT_CA1915AA_AS',
    'ENT_CA1916AA_AS', 'ENT_CA1917AA_AS', 'ENT_CA1918AA_AS', 'ENT_CA1919AA_AS'
]

def should_exclude(seq_id):
    for pattern in exclude_patterns:
        if pattern in seq_id:
            return True
    return False

def calculate_conservation(alignment_file):
    sequences = list(SeqIO.parse(alignment_file, 'fasta'))
    if not sequences:
        return None, 0, 0

    aln_length = len(sequences[0].seq)
    n_seqs = len(sequences)

    entropies = []
    for pos in range(aln_length):
        column = [str(seq.seq[pos]).upper() for seq in sequences]
        counts = Counter(column)
        total = sum(counts.values())

        entropy = 0
        for count in counts.values():
            if count > 0:
                p = count / total
                entropy -= p * np.log2(p)

        max_entropy = np.log2(min(4, total))
        norm_entropy = entropy / max_entropy if max_entropy > 0 else 0
        entropies.append(norm_entropy)

    return 1 - np.mean(entropies), n_seqs, aln_length

# Process DNA alignments
msa_dir = 'output/msa/dna_alignments'
results = []

for msa_file in os.listdir(msa_dir):
    if not msa_file.endswith('.fasta'):
        continue

    filepath = os.path.join(msa_dir, msa_file)
    gene_name = msa_file.replace('_dna.fasta', '').replace('.fasta', '')

    sequences = list(SeqIO.parse(filepath, 'fasta'))
    original = len(sequences)

    filtered = [seq for seq in sequences if not should_exclude(seq.id)]
    removed = original - len(filtered)

    if removed > 0:
        SeqIO.write(filtered, filepath, 'fasta')
        print(f"  {gene_name}: removed {removed} sequences")

    conservation, n_seqs, aln_len = calculate_conservation(filepath)
    results.append({
        'gene': gene_name,
        'n_sequences': n_seqs,
        'alignment_length': aln_len,
        'conservation_score': conservation
    })

df = pd.DataFrame(results)
df.to_csv('output/operon_conservation_metrics.csv', index=False)

# Update manuscript_stats.txt
with open('manuscript_stats.txt', 'w') as f:
    f.write(f"Operon Assembly Extraction Statistics\n")
    f.write(f"=====================================\n")
    f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    f.write(f"Genomes analyzed: 8,596\n")
    f.write(f"Operon genes extracted: {len(df)}\n")
    f.write(f"Average sequences per gene: {df['n_sequences'].mean():.0f}\n\n")
    f.write(f"Per-gene statistics:\n")
    for _, row in df.iterrows():
        f.write(f"  {row['gene']}: {row['n_sequences']} sequences, conservation={row['conservation_score']:.3f}\n")

print(f"\nUpdated operon conservation metrics")
print(df.to_string(index=False))
PYTHON
```

---

### Phase 2: Re-run Downstream Analyses

These steps need to be re-run because they compute new metrics from the filtered inputs:

#### Step 2.1: Diversity Analysis (06_diversity_analysis)

```bash
cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/06_diversity_analysis
conda activate efs_diversity

# Re-run comparative analysis (reads from 04 and 05 outputs)
python comparative_analysis.py

# Or if there's a manuscript_numbers.py:
python manuscript_numbers.py > manuscript_stats.txt 2>/dev/null || echo "Using comparative_analysis.py output"
```

---

#### Step 2.2: dN/dS Analysis (07_dnds_analysis) - LONG JOB

**Note**: This is computationally intensive. Consider if needed.

```bash
cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/07_dnds_analysis

# Submit as SLURM job (uses filtered MSAs from steps 04/05)
sbatch run_dnds_analysis.sh
```

---

#### Step 2.3: Codon Variation (07b) - After dN/dS completes

```bash
cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/07b_codon_variation_analysis
sbatch run_codon_variation.sh
```

---

#### Step 2.4: Variation Visualization (07c)

```bash
cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/07c_visualize_variation_analysis
sbatch run_visualization.sh --genes ptsA ptsB ptsC ptsD frpC glpC fruR
```

---

#### Step 2.5: Start Site Analysis (08_start_site_analysis)

**Action**: Re-run with filtering

```bash
cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/08_start_site_analysis
conda activate efs_diversity

# The script reads from prokka outputs - need to filter
# Option 1: Modify the script to use valid_as_barcodes.txt
# Option 2: Run and then filter output

sbatch run_start_site_analysis.sh

# After completion, filter if needed:
python3 << 'PYTHON'
import pandas as pd

with open('../valid_as_barcodes.txt') as f:
    valid_ids = set(line.strip() for line in f)

# Filter start_site_summary.tsv if it exists
try:
    df = pd.read_csv('output/start_site_summary.tsv', sep='\t')
    df['genome_id'] = df['genome'].str.replace('.result', '').str.replace('.genome', '')
    df_filtered = df[df['genome_id'].isin(valid_ids)]
    df_filtered.drop(columns=['genome_id']).to_csv('output/start_site_summary.tsv', sep='\t', index=False)
    print(f"Filtered: {len(df)} -> {len(df_filtered)} rows")
except Exception as e:
    print(f"Note: {e}")
PYTHON
```

---

#### Step 2.6: Operon Order Analysis (09_operon_order_analysis)

```bash
cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/09_operon_order_analysis

# Re-run (reads from BLAST results which are now filtered)
sbatch run_pipeline.sh

# Generate updated stats
python manuscript_numbers.py > output/manuscript_stats.txt 2>/dev/null || true
```

---

#### Step 2.7: Full Operon Mapping (13e_full_operon_mapping)

```bash
cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/13e_full_operon_mapping

# This maps to assemblies - need to filter the assembly list
# Create filtered assembly list
python3 << 'PYTHON'
import os
with open('../valid_as_barcodes.txt') as f:
    valid_ids = set(line.strip() for line in f)

# The script likely scans Efs_assemblies - need to check how it works
# May need to modify run_full_mapping.sh to use a filtered list
print(f"Valid assemblies: {len(valid_ids)}")
print("Check run_full_mapping.sh to see how it selects assemblies")
PYTHON

# Then run the mapping
sbatch run_full_mapping.sh
# ... follow with other steps per original PLAN
```

---

## Execution Order Summary

```
Phase 1 (Filter only - fast):
  1.1 Prokka stats
  1.2 BLAST results
  1.3 Core gene MSAs (parallel possible)
  1.4 Operon MSAs (parallel possible)

Phase 2 (Re-run - variable time):
  2.1 Diversity analysis (fast, ~minutes)
  2.2 dN/dS analysis (slow, ~hours-days)
  2.3 Codon variation (after 2.2)
  2.4 Variation viz (fast)
  2.5 Start site analysis (fast)
  2.6 Operon order (fast)
  2.7 Full mapping (medium)
```

---

## Quick Start: Run Phase 1

```bash
cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis
conda activate efs_diversity

# Run all Phase 1 steps (copy-paste each section from above)
# Or create a script:

cat > run_phase1_filtering.sh << 'SCRIPT'
#!/bin/bash
set -e

echo "Phase 1: Filtering existing outputs..."
echo ""

# Step 1.1
echo "=== Step 1.1: Prokka stats ==="
cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/01_prokka_annotation
# ... (python code from above)

# Continue with other steps...
SCRIPT

chmod +x run_phase1_filtering.sh
```

---

## Verification After Completion

```bash
# Check all manuscript_stats.txt show 8,596
grep -rh "8,596\|8596" */manuscript_stats.txt */output/manuscript_stats.txt 2>/dev/null

# Verify no excluded IDs remain in MSAs
for pattern in ENT_CA1913AA_AS ENT_CA1914AA_AS ENT_CA1915AA_AS ENT_CA1916AA_AS ENT_CA1917AA_AS ENT_CA1918AA_AS ENT_CA1919AA_AS; do
    echo -n "$pattern: "
    grep -r "$pattern" 04_core_gene_analysis/output/core_gene_alignments/ 05_operon_assembly_extraction/output/msa/ 2>/dev/null | wc -l
done
```

---

## Notes

1. **MSA filtering preserves alignment**: Removing sequences from an aligned FASTA doesn't break the alignment - remaining sequences stay aligned.

2. **Conservation scores will change slightly**: With 7 fewer sequences (0.08%), scores may shift by small amounts.

3. **dN/dS is the slowest step**: Consider if it's worth re-running. The 7 excluded sequences likely have minimal impact.

4. **Backup before filtering**:
   ```bash
   cp -r 04_core_gene_analysis/output/core_gene_alignments{,_backup_8603}
   cp -r 05_operon_assembly_extraction/output/msa{,_backup_8603}
   ```
