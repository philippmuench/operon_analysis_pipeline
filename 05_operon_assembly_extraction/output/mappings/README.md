# Operon Sequence Extraction Strategies

This directory contains results from 4 different approaches to extract operon gene sequences from bacterial genomes.

## Strategy Overview

### **Strategy A**: `aa_nt_mapping/prokka/` 
**Protein → Nucleotide mapping using Prokka annotations**
- **Query**: Reference protein sequences (7 operon genes)
- **Target**: Prokka-annotated genomes (.fna files)
- **BLAST**: `tblastn` (protein vs nucleotide)
- **Purpose**: Find genes using protein similarity, extract nucleotide sequences
- **Use case**: Gene boundary analysis and validation

### **Strategy B**: `nt_nt_mapping/prokka_genome/`
**Nucleotide → Nucleotide mapping using Prokka genomes**
- **Query**: Reference nucleotide sequences (7 operon genes)  
- **Target**: Prokka-annotated genomes (.fna files)
- **BLAST**: `blastn` (nucleotide vs nucleotide)
- **Purpose**: Direct nucleotide sequence matching
- **Use case**: Standard homology search

### **Strategy C**: `nt_nt_mapping/prokka_variants/`
**Reverse mapping to capture sequence variants**
- **Query**: Prokka-predicted genes (.ffn files from each genome)
- **Target**: Reference operon gene database
- **BLAST**: `blastn` (predicted genes vs reference)
- **Purpose**: "Which predicted genes are operon homologs?" → captures natural variants
- **Use case**: Phylogenetic analysis (preserves sequence diversity)

### **Strategy D**: `aa_nt_mapping/assemblies/`
**Direct assembly mapping (no annotation dependency)**
- **Query**: Reference protein sequences (7 operon genes)
- **Target**: Raw genome assemblies (.fasta.gz files)
- **BLAST**: `tblastn` (protein vs nucleotide)
- **Purpose**: Find genes directly in assemblies, bypass annotation issues
- **Use case**: Primary strategy for downstream phylogenetic analysis

## Key Differences

| Strategy | What's being mapped | Direction | Result |
|----------|-------------------|-----------|---------|
| A | Reference proteins → Prokka genomes | Forward | Protein-guided gene finding |
| B | Reference genes → Prokka genomes | Forward | Direct sequence matching |
| C | Genome genes → Reference database | **Reverse** | Natural sequence variants |
| D | Reference proteins → Raw assemblies | Forward | Annotation-independent |

## Directory Contents

Each strategy folder contains:
- `sequences/` - Extracted FASTA files
- `msa/` - Multiple sequence alignments  
- `plots/` - Conservation plots (SNP-based)
- `enhanced_plots/` - Advanced conservation plots (Shannon entropy + sequence logos)

## Recommended Usage

- **Strategy A**: Gene boundary validation
- **Strategy D**: Primary downstream analysis (phylogenetics, evolution)
- **Strategy C**: Sequence diversity analysis
- **Strategy B**: Standard comparative genomics
