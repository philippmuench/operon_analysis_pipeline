# Step 1: Prokka Annotation

This directory contains scripts for annotating E. faecalis genomes using Prokka.

## Input
- **Source**: `../Efs_assemblies/` - 8,587 compressed E. faecalis genome assemblies

## Output
- **Production mode**: `../prokka_output/` - One subdirectory per genome
- **Test mode**: `output/prokka_results/` - Results for first 50 genomes only
- **Each genome directory contains**:
  - `.gff` - Gene annotations
  - `.faa` - Protein sequences (amino acids)
  - `.ffn` - Nucleotide sequences (genes)
  - `.fna` - Nucleotide sequences (contigs)
  - `.gbk` - GenBank format
  - `.tsv` - Feature table
  - `.txt` - Statistics
  - `.log` - Prokka log file

## Usage
```bash
# Run full analysis on all 8,587 genomes
sbatch run_prokka.sh

# Test mode: process only first 50 genomes
sbatch run_prokka.sh --test
```

This SLURM array job processes genomes in batches of 100, running up to 20 batches in parallel. The `--test` parameter is useful for testing the pipeline on a smaller subset of genomes.

## Monitoring
```bash
# Check progress
./check_progress.sh

# Validate all outputs are complete and non-empty
./validate_prokka_outputs.sh

# Rerun any failed/incomplete genomes
sbatch rerun_failed_genomes.sh [incomplete_genomes_list.txt]
```

## Validation

The `validate_prokka_outputs.sh` script performs comprehensive checks:
- Verifies all 8 expected output files exist for each genome
- Checks that all files are non-empty
- Reports detailed statistics and identifies problematic genomes
- Generates a list of incomplete genomes for reprocessing
- Creates a timestamped validation report

## Troubleshooting

If genomes fail or have incomplete outputs:
1. Run `./validate_prokka_outputs.sh` to identify issues
2. Review the generated log file for specific problems
3. Use `sbatch rerun_failed_genomes.sh` to reprocess failed genomes
4. The rerun script will automatically validate outputs when complete