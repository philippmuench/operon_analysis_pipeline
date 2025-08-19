### Start-site analysis for operon genes

This module investigates start-codon choice discrepancies between reference operon genes (`02_reference_operon_extraction/output/operon_genes_nt.fasta`) and Prokka/Prodigal calls across assemblies. It scans each called CDS for upstream, in-frame alternative start codons (ATG/GTG/TTG), checks for Shine–Dalgarno-like RBS motifs at appropriate spacings, and summarizes why a shorter downstream start may have been chosen.

#### Inputs
- `--prokka_dir` directory with per-genome Prokka outputs (each subdir contains `.gff` and `.fna`).
- `--gene_reference_fasta` multi-FASTA of operon reference genes (default points to `02_reference_operon_extraction/output/operon_genes_nt.fasta`). Used to determine the set of gene names to analyze.
- Optional: `--genome_list` file with one genome id per line to limit the analysis.

#### Outputs (in `08_start_site_analysis/output` by default)
- `start_site_summary.tsv`: One row per genome × gene with fields including:
  - genome_id, gene, contig, strand, start, end
  - start_codon (from sequence), `start_type`/`partial`/`rbs_motif`/`rbs_spacer` (from GFF when present)
  - upstream_candidate_found, upstream_codon, upstream_offset_nt
  - upstream_has_rbs, upstream_rbs_seq, upstream_rbs_spacer_nt, upstream_rbs_mismatches
  - classification (e.g., Alt_upstream_with_RBS, Upstream_no_RBS, No_upstream, Partial_gene, Same_start)
- `contexts/` (optional): per-case sequence context snippets around the chosen and alternative starts.
- `report.html` (optional): compact HTML with highlighted starts/RBS for a sample of interesting cases.

#### Heuristics
- Scan up to `--max_upstream` nt upstream (default 120 nt) in the gene’s frame for ATG/GTG/TTG that extend the ORF without encountering an in-frame stop (TAA/TAG/TGA).
- RBS-like motifs checked within `--rbs_min_spacer`..`--rbs_max_spacer` (default 4..13 nt) upstream of each candidate start. Motifs tested: `AGGAGG`, `GGAGG`, `AGGA`, `GGAG`, `GAGG`. Up to 1 mismatch is allowed.

#### Usage
```bash
cd 08_start_site_analysis

# Run locally on a subset
python analyze_start_sites.py \
  --prokka_dir ../prokka_output \
  --gene_reference_fasta ../02_reference_operon_extraction/output/operon_genes_nt.fasta \
  --output_dir ./output \
  --max_workers 8

# Or submit to SLURM for the full dataset
sbatch run_start_site_analysis.sh
```

#### Notes
- The script is robust to missing attributes; if `start_type`/`rbs_motif` are absent in GFF, they are inferred from sequence where possible.
- For minus-strand genes, sequences are analyzed in the coding orientation (reverse-complement and reindexed so position 0 is the called start codon).
- This module does not change calls; it explains likely reasons why upstream starts were not chosen.


