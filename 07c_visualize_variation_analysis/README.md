# Codon Variation Visualisation

## Overview

This step packages a lightweight workflow that turns operon multiple-sequence
alignments into ready-to-plot summaries of codon-level polymorphisms. It uses
the same counting logic as the `07b` selection analysis so the numbers line up
exactly, but adds a compact per-codon table and a colour-coded plot that can be
used directly in figures.

## Why this step?

- Highlight which codon positions contribute to the `codon_variable_sites` count.
- Distinguish synonymous-only, non-synonymous-only, and mixed sites.
- Provide quick context (allele counts, amino-acid states) for figure captions.

## Inputs

- Codon-aligned FASTA files from `05_operon_assembly_extraction/output/msa/dna_alignments`.
- By default, the script processes `ptsA_aligned.fasta`, but any gene name may be
  supplied (e.g. `--genes ptsA fruR`).

## Usage

### Batch mode (recommended on the cluster)

```bash
sbatch run_visualization.sh --genes ptsA
```

Use the `--genes` flag to supply one or more gene names; they must match the
`*_aligned.fasta` files in the alignment directory.

### Manual run

```bash
conda activate efs_diversity
python visualize_codon_variation.py \
    --genes ptsA \
    --alignment-dir ../05_operon_assembly_extraction/output/msa/dna_alignments \
    --output-dir output \
    --report \
    --focus-codons 145 \
    --max-examples 3 \
    --window-size 10
```

## Outputs

All results are written beneath `output/`:

- `tables/<gene>_variable_codons.csv` – per-codon summary table with categories,
  allele counts, and amino-acid states.
- `plots/<gene>_variable_codons.pdf` – scatter plot marking each variable codon,
  coloured by synonymous/non-synonymous status.
- `plots/focus_sites/<gene>_codon_<index>_profile.pdf` – bar chart of allele
  counts for highlighted codon positions.
- `plots/focus_sites/examples/<gene>_codon_<index>_examples.pdf` – raw codon
  snapshots showing example sequences (colour-coded by nucleotide and amino acid).
- `plots/focus_sites/windows/<gene>_codon_<index>_window.pdf` – heatmap of amino
  acid usage across a configurable window (±`--window-size` codons) with per-site
  codon diversity summaries.
- `codon_variation_highlights.txt` (optional, `--report`) – text summary for quick
  reference or manuscript snippets.

## Notes

- Codons containing gaps or ambiguous bases are skipped, so indels do not create
  spurious variation.
- The colours are consistent with the `07b` plots (green = synonymous-only,
  red = non-synonymous-only, purple = mixed).
- Point sizes scale with the number of distinct codon alleles observed; larger
  markers flag richer diversity at that position.
