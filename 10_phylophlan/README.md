# PhyloPhlan Analysis for Enterococcus faecalis

This directory contains the PhyloPhlan pipeline for phylogenetic analysis of *Enterococcus faecalis* isolates.

## Overview

PhyloPhlan is used to build high-resolution phylogenetic trees based on core gene alignments. This analysis processes a subset of 100 E. faecalis genomes for testing purposes.

## Prerequisites

- PhyloPhlan 3.0+
- Diamond
- MAFFT
- TrimAl
- FastTree
- RAxML
- GraphLan (optional, for visualization)

## Setup Steps

### Step 1: Request Interactive Resources

Request an interactive SLURM session with sufficient resources:

```bash
srun --cpus-per-task=4 --mem=10G --time=04:00:00 --pty bash
```

### Step 2: Setup PhyloPhlan Database

Download and setup the species-specific database for *Enterococcus faecalis*:

```bash
phylophlan_setup_database \
    -g s__Enterococcus_faecalis \
    -o . \
    --verbose 2>&1 | tee logs/phylophlan_setup_database.log
```

This command:
- Downloads the core gene markers for *E. faecalis*
- Creates a Diamond database for sequence mapping
- Outputs to the current directory

### Step 3: Generate Configuration File

Create the configuration file that specifies the tools and parameters for the analysis:

```bash
phylophlan_write_config_file \
    -o isolates_config.cfg \
    -d a \
    --force_nucleotides \
    --db_aa diamond \
    --map_aa diamond \
    --map_dna diamond \
    --msa mafft \
    --trim trimal \
    --tree1 fasttree \
    --tree2 raxml
```

Configuration parameters:
- `-d a`: Analyze amino acid sequences
- `--force_nucleotides`: Force nucleotide-based analysis
- `--db_aa diamond`: Use Diamond for protein database searches
- `--msa mafft`: Use MAFFT for multiple sequence alignment
- `--trim trimal`: Use TrimAl for alignment trimming
- `--tree1 fasttree`: Use FastTree for initial tree construction
- `--tree2 raxml`: Use RAxML for refined tree construction

### Step 4: Run PhyloPhlan Analysis

Execute the main phylogenetic analysis pipeline:

```bash
phylophlan \
    -i input_isolates \
    -o output_isolates \
    -d s__Enterococcus_faecalis \
    --trim greedy \
    --not_variant_threshold 0.99 \
    --remove_fragmentary_entries \
    --fragmentary_threshold 0.67 \
    --min_num_entries 50 \
    -t a \
    -f isolates_config.cfg \
    --diversity low \
    --force_nucleotides \
    --nproc 4 \
    --verbose 2>&1 | tee logs/phylophlan__output_isolates.log
```

Key parameters:
- `-i input_isolates`: Input directory containing genome sequences
- `-o output_isolates`: Output directory for results
- `--trim greedy`: Use greedy trimming strategy
- `--not_variant_threshold 0.99`: Threshold for identifying variants (99% identity)
- `--remove_fragmentary_entries`: Remove incomplete sequences
- `--fragmentary_threshold 0.67`: Consider sequences <67% complete as fragmentary
- `--min_num_entries 50`: Minimum number of genomes required for marker inclusion
- `--diversity low`: Optimized for low diversity datasets (same species)
- `--nproc 4`: Use 4 CPU cores

## Running the Pipeline

For batch processing, use the provided script:

```bash
sbatch run_phylophlan.sh
```

This script submits a SLURM job that runs the complete PhyloPhlan analysis.

## Output Files

The analysis produces several key outputs in the `output_isolates/` directory:

- `input_isolates.tre`: Initial phylogenetic tree (FastTree)
- `input_isolates_resolved.tre`: Resolved tree with polytomies handled
- `RAxML_bestTree.input_isolates_refined.tre`: Best tree from RAxML analysis
- `RAxML_result.input_isolates_refined.tre`: Final refined phylogenetic tree
- `RAxML_info.input_isolates_refined.tre`: RAxML run information
- `RAxML_log.input_isolates_refined.tre`: RAxML execution log
- `input_isolates_concatenated.aln`: Concatenated multiple sequence alignment
- `input_isolates_concatenated.aln.reduced`: Reduced alignment after trimming

## Optional: Tree Visualization with GraphLan

To create annotated circular phylogenetic tree visualizations:

```bash
graphlan_annotate.py \
    --annot graphlan/isolates_annotation.txt \
    output_isolates/RAxML_bestTree.input_isolates.tre \
    graphlan/isolates_annotated.xml
```

## Directory Structure

```
10_phylophlan/
├── README.md                       # This file
├── run_phylophlan.sh              # Main execution script
├── isolates_config.cfg            # PhyloPhlan configuration
├── input_isolates/                # Input genome sequences (100 test genomes)
├── output_isolates/               # Analysis outputs
│   ├── *.tre                      # Phylogenetic trees
│   └── tmp/                       # Temporary processing files
├── s__Enterococcus_faecalis/      # Species database files
└── logs/                          # Execution logs
```

## Notes

- Currently configured for 100 test genomes
- Full analysis would process all available E. faecalis genomes
- Tree files are in Newick format, compatible with most phylogenetic visualization tools
- The analysis focuses on core gene conservation across isolates