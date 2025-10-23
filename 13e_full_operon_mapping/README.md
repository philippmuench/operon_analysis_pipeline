# Step 13e: Full Operon Mapping

Compare how the complete operon sequence (from GenBank ORIGIN) maps to each
assembly for both the legacy reference and the SNP-aware update.

## Usage
```bash
python compare_full_operon_mapping.py \
    --legacy-gb ../02_reference_operon_extraction/operon.gb \
    --updated-gb ../13a_new_reference_operon_extraction/FL_operon_with_SNPs.gb \
    --assemblies-dir ../../Efs_assemblies \
    --output output_full_mapping --threads 4
```

Add `--max-assemblies N` for a quick dry run. Results include:
- `full_operon_mapping.tsv` – per-assembly hit metrics (identity/coverage)
- `summary.txt` – aggregated hit counts (total vs high-quality hits)

Adjust `--min-coverage`/`--min-identity` thresholds as needed.
