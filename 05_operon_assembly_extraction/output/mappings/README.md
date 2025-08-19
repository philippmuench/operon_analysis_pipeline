# Operon Extraction Strategies

Four different approaches to extract operon gene sequences:

## Strategies

| Strategy | Query → Target | Method | Purpose |
|----------|----------------|---------|---------|
| **A** `aa_nt_mapping/prokka/` | Reference proteins → Prokka genomes | tblastn | Gene boundary analysis |
| **B** `nt_nt_mapping/prokka_genome/` | Reference genes → Prokka genomes | blastn | Direct sequence matching |
| **C** `nt_nt_mapping/prokka_variants/` | Genome genes → Reference database | blastn (reverse) | Capture sequence variants |
| **D** `aa_nt_mapping/assemblies/` | Reference proteins → Raw assemblies | tblastn | Primary analysis (no annotation) |

## Key Insight - Strategy C
Strategy C does **reverse BLAST**: asks "which genome-predicted genes match our reference?" to capture natural sequence variants.

## Folder Contents
- `sequences/` - Extracted sequences
- `msa/` - Alignments  
- `plots/` - Conservation plots
- `enhanced_plots/` - Shannon entropy + sequence logos
