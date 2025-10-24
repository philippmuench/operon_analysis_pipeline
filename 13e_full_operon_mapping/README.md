# Step 13e: Full Operon Mapping

Compare how the complete operon sequence (from GenBank ORIGIN) maps to each
assembly for both the legacy reference and the SNP-aware update.

## Usage
```bash
python compare_full_operon_mapping.py \
    --legacy-gb ../02_reference_operon_extraction/operon.gb \
    --updated-gb ../13a_new_reference_operon_extraction/FL_operon_with_SNPs.gb \
    --assemblies-dir ../../Efs_assemblies \
    --output output_full_mapping --threads 4 --min-coverage 80 --min-identity 90 --save-raw
```

Add `--max-assemblies N` for a quick dry run. Results include:
- `full_operon_mapping.tsv` – per-assembly hit metrics (identity/coverage and aligned strings)
- `summary.txt` – aggregated hit counts (total vs high-quality hits)
- `raw_blast/legacy|updated/*.tsv` – raw BLAST tabular output (when `--save-raw` is set) for downstream SNP/gap analyses

Adjust `--min-coverage`/`--min-identity` thresholds as needed.

## Post-processing raw BLAST output

Once the mapping job completes, convert the raw BLAST alignments into per-nucleotide
mismatch/deletion summaries via:

```bash
sbatch run_process_raw.sh
```

Outputs are written to `output_full_run/position_stats/`, including
`legacy_position_summary.tsv`, `updated_position_summary.tsv`, and a combined table
that can be used to visualise conservation or gap rates across the operon.

## Visualising position-level metrics

After generating the position summaries, create operon-wide plots via:

```bash
sbatch run_position_plots.sh
```

Plots (mismatch rate, deletion rate, coverage, and their deltas) are saved in
`output_full_run/position_stats/plots/`.

## Orientation summary

Run the following to tabulate alignment direction (query forward vs reverse) for each dataset:

```bash
python - <<'PY'
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

raw_base = Path('output_full_run/raw_blast')
results = []
for label in ['legacy','updated']:
    files = sorted((raw_base / label).glob('*.tsv'))
    forward = reverse = total = 0
    for file in files:
        best_line = None
        best_bitscore = None
        for line in file.read_text().splitlines():
            if not line:
                continue
            parts = line.split('\t')
            bitscore = float(parts[11])
            if best_bitscore is None or bitscore > best_bitscore:
                best_bitscore = bitscore
                best_line = parts
        if best_line is None:
            continue
        qstart = int(best_line[6])
        qend = int(best_line[7])
        total += 1
        if qend >= qstart:
            forward += 1
        else:
            reverse += 1
    results.append({'dataset': label, 'forward': forward, 'reverse': reverse, 'total': total})

summary = pd.DataFrame(results)
summary_path = Path('output_full_run/orientation_breakdown.tsv')
summary.to_csv(summary_path, sep='\t', index=False)
print(summary)

summary.set_index('dataset')[['forward','reverse']].plot(kind='bar', stacked=True, color=['#4C72B0','#C44E52'])
plt.ylabel('Assemblies (best hit)')
plt.title('Full operon alignment orientation')
plt.tight_layout()
plt.savefig('output_full_run/orientation_breakdown.png')
plt.close()
PY
```

This produces `orientation_breakdown.tsv` and `orientation_breakdown.png` summarising forward vs reverse alignments.

Zoomed deletion-rate plots (highlighting regions where deletion rate exceeds 0.015) are automatically generated as `*_deletion_rate_zoom.png`. Adjust the threshold via `--deletion-threshold` if needed.
