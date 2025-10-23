# Step 13b: BLAST Search With SNP-Aware Reference

This directory mirrors the Stage 03 BLAST pipeline but points to the SNP-aware operon
reference generated in `../13a_new_reference_operon_extraction/output/`.

## Input
- Query proteins: `../13a_new_reference_operon_extraction/output/operon_genes_protein.fasta`
- Query coding (nt): `../13a_new_reference_operon_extraction/output/operon_genes_nt.fasta`
- Query regulatory: `../13a_new_reference_operon_extraction/output/operon_regulatory_nt.fasta`
- Genome databases: `../../Efs_assemblies/*.fasta.gz`
- Prokka annotations: `../01_prokka_annotation/output/prokka_results/*/`

## Usage

### Run BLAST search batches (array job)
```bash
sbatch --array=1-86 run_pipeline.sh search
```

### Process and summarise results once searches finish
```bash
sbatch run_pipeline.sh all
```

### Direct execution (example)
```bash
python blast_pipeline.py --mode all --threads 60 \
    --ref-dir ../13a_new_reference_operon_extraction/output
```

Outputs are written to the local `output/` directory, mirroring the Stage 03 layout.
