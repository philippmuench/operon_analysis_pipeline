# Step 2: Operon Gene Extraction

This directory contains scripts for extracting operon genes from the reference GenBank file.

## Input
- **Source**: `operon.gb` - GenBank file containing the fructoselysine/glucoselysine operon (local file)

## Output
- `output/operon_genes.tsv` - Table with gene information (locus tags, names, products, positions)
- `output/operon_genes.fasta` - Protein sequences for BLAST searches

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