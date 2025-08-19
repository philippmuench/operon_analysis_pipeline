# Step 5: Operon Assembly Extraction - Restart Guide

## Quick Restart Commands

### 1. Full Pipeline Restart (Recommended)
```bash
# Clean restart from the beginning
cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/05_operon_assembly_extraction

# Submit via SLURM (recommended for full pipeline)
sbatch run_operon_extraction.sh

# Monitor progress
squeue -u $USER
tail -f operon_extraction_*.out
```

### 2. Manual Step-by-Step Restart
```bash
cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/05_operon_assembly_extraction

# Step 1: Extract sequences using BLAST coordinates from step 3
python extract_operon_sequences.py \
    --prokka_dir ../01_prokka_annotation/output/prokka_results \
    --blast_dir ../03_blast_search/output/blast_results \
    --output_dir output/sequences \
    --min_identity 90 \
    --min_coverage 80

# Step 2: Create MSAs
python create_msa.py \
    --coding-sequences output/sequences \
    --output-dir output/msa \
    --threads 8

# Step 3: Extract promoter sequences
python extract_noncoding_sequences.py \
    --blast-dir ../03_blast_search/output/blast_results \
    --output-dir output/noncoding_sequences \
    --min-identity 70

# Step 4: Create conservation plots
python create_gene_conservation_plots.py \
    --msa-dir output/msa/dna_alignments \
    --output-dir output/plots/gene_conservation

python create_promoter_plot_with_pribnow.py \
    --blast-based \
    --output-dir output/plots
```

### 3. Start from Specific Step
```bash
# Start from step 2 (MSA creation) if sequences already extracted
sbatch run_operon_extraction.sh --start-step 2

# Start from step 3 (promoter analysis)
sbatch run_operon_extraction.sh --start-step 3

# Start from step 4 (plotting only)
sbatch run_operon_extraction.sh --start-step 4
```

## Prerequisites

### Required Input from Previous Steps:
- ✅ **Step 3 BLAST results**: `../03_blast_search/output/blast_results/`
- ✅ **Prokka annotations**: `../01_prokka_annotation/output/prokka_results/`
- ✅ **Reference sequences**: `../02_reference_operon_extraction/output/`

### Check Prerequisites:
```bash
# Verify BLAST results exist
ls ../03_blast_search/output/blast_results/*.txt | wc -l
# Should show thousands of BLAST result files

# Verify Prokka annotations exist
ls ../01_prokka_annotation/output/prokka_results/ | wc -l
# Should show ~8,587 genome directories

# Verify reference sequences exist
ls ../02_reference_operon_extraction/output/operon_genes_*.fasta
# Should show reference FASTA files
```

## Resource Requirements

### SLURM Resources (recommended):
```bash
sbatch -c 16 --mem=24G --time=6:00:00 run_operon_extraction.sh
```

### For manual runs:
- **CPU**: 8-16 cores for MAFFT alignments
- **Memory**: 16-24GB for large MSAs
- **Time**: 2-4 hours depending on data size
- **Storage**: High-performance local scratch (`/vol/tmp`) for MAFFT

## Expected Outputs

After successful completion, you should have:

```
output/
├── sequences/                   # ~7 gene FASTA files
├── msa/
│   ├── dna_alignments/         # DNA MSAs for each gene
│   ├── protein_alignments/     # Protein MSAs (if generated)
│   └── noncoding_alignments/   # Promoter MSAs
├── plots/                      # Conservation plots
├── noncoding_sequences/        # Promoter sequences
└── extraction_pipeline_summary.txt
```

## Troubleshooting

### Common Issues:

1. **"BLAST results not found"**
   ```bash
   # Check if step 3 completed successfully
   ls ../03_blast_search/output/blast_results/ | head -10
   ```

2. **"No sequences extracted"**
   ```bash
   # Check BLAST result quality
   head ../03_blast_search/output/blast_results/ENT_*.txt
   # Verify identity/coverage thresholds aren't too strict
   ```

3. **"MAFFT alignment failed"**
   ```bash
   # Check available memory and CPU
   # Ensure /vol/tmp is available and writable
   ls -la /vol/tmp/
   ```

4. **"Permission denied on /vol/tmp"**
   ```bash
   # Create user-specific temp directory
   export MAFFT_TMPDIR="/tmp/mafft_$$"
   mkdir -p $MAFFT_TMPDIR
   ```

## Clean Restart (if needed)

```bash
# Remove all outputs to start fresh
rm -rf output/*
mkdir -p output/{sequences,msa,noncoding_sequences,plots}

# Then restart with full pipeline
sbatch run_operon_extraction.sh
```

## Dependencies

This step depends on **COMPLETED** outputs from:
- ✅ Step 1: Prokka annotation
- ✅ Step 2: Reference extraction  
- ✅ Step 3: BLAST search

**Note**: This step does NOT run new BLAST searches - it uses existing results from step 3!
