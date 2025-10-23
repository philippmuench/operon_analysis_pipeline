# Step 13d: Operon Extraction Comparison

Compare the `extraction_pipeline_summary.txt` outputs from two pipeline runs (e.g.
legacy Step 05 vs. the SNP-aware Step 13c).

## Usage
```
python compare_extraction_summaries.py \
    ../05_operon_assembly_extraction/output/extraction_pipeline_summary.txt \
    ../13c_new_operon_assembly_extraction/output/extraction_pipeline_summary.txt \
    --labels legacy updated
```

Outputs are written to `output/`:

- `sequences_comparison.tsv` – per-gene sequence counts with differences
- `alignments_comparison.tsv` – per-gene alignment lengths and sequence counts
- `top_conservation_comparison.tsv` – optional table if top scores are listed
- `summary.txt` – quick description with input paths

### Conservation metrics comparison
```
python compare_conservation_metrics.py --legacy-label legacy --updated-label updated
```

This writes tables under `output/metrics_comparison/` and scatter plots in
`output/metrics_comparison/plots/` showing updated (Step 13) metrics versus legacy (Step 05).
