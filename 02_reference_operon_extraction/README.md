# Step 2: Reference Operon Extraction (from GenBank)

This directory contains scripts for parsing the reference operon from a curated GenBank file and exporting gene and promoter sequences for downstream BLAST queries.

## Input
- **Source**: `operon.gb` - GenBank file containing the fructoselysine/glucoselysine operon (local file)

## Output
- `output/operon_genes.tsv` - Table with gene information (locus tags, names, products, positions)
- `output/operon_genes_nt.fasta` - Nucleotide sequences for CDS
- `output/operon_genes_protein.fasta` - Protein sequences for BLAST searches
- `output/operon_noncoding_nt.fasta` - Promoter and other non-coding sequences

## Usage
```bash
python extract_operon_genes.py
```

## Operon Genes
The fructoselysine/glucoselysine operon contains 7 genes:
1. Fructoselysine-6-phosphate deglycase
2. Glucoselysine-6-phosphate deglycase  
3. PTS system EIID component
4. PTS system EIIC component
5. PTS system EIIB component
6. PTS system EIIA component
7. Sigma-54 dependent transcriptional regulator

## Manuscript Statistics
Generate statistics for the manuscript after extracting reference sequences:
```bash
# Generate with SLURM (recommended)
sbatch run_manuscript_stats.sh

# Or run directly (outputs to console)
python manuscript_numbers.py

# Save statistics to file
python manuscript_numbers.py manuscript_stats.txt
```

The statistics include:
- Reference genome source information
- Number and length of extracted protein sequences
- Number and length of extracted nucleotide sequences
- Operon gene composition details

## Notes
- These exported sequences are used by Step 03 (`../03_blast_search/`) as BLAST queries.
- If the GenBank content changes, re-run this step before repeating the BLAST searches.