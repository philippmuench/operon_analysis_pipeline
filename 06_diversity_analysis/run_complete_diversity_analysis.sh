#!/bin/bash
#SBATCH --job-name=diversity_complete
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=output/diversity_complete_%j.out
#SBATCH --error=output/diversity_complete_%j.err

################################################################################
# COMPLETE DIVERSITY ANALYSIS PIPELINE
# 
# This script runs the full diversity analysis including:
# 1. Coding gene sequence extraction and MSA
# 2. Non-coding sequence (promoter) extraction and MSA  
# 3. Conservation analysis for both sequence types
# 4. Gap analysis across all alignments
# 5. All visualization plots
#
# All outputs are saved to output/diversity_analysis/
################################################################################

echo "=========================================="
echo "COMPLETE DIVERSITY ANALYSIS PIPELINE"
echo "=========================================="
echo "Starting at $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $SLURMD_NODENAME"

# Create output directory
mkdir -p output

# Load conda environment
eval "$(/home/pmuench/miniconda3/bin/conda shell.bash hook)"
conda activate efs_diversity

# Check if conda environment loaded successfully
if [ $? -ne 0 ]; then
    echo "Error: Failed to activate conda environment"
    exit 1
fi

echo "Using Python: $(which python)"
echo "Using MAFFT: $(which mafft)"

# Set number of threads
THREADS=${SLURM_CPUS_PER_TASK:-4}
echo "Using $THREADS threads"

################################################################################
# STEP 1: EXTRACT CODING GENE SEQUENCES (FIXED PROKKA PATHS)
################################################################################
echo ""
echo "STEP 1: Extracting coding gene sequences..."
echo "----------------------------------------"
python extract_operon_sequences.py \
    --prokka_dir ../01_prokka_annotation/output/prokka_results \
    --blast_dir ../03_blast_search/output/blast_results \
    --output_dir output/operon_sequences \
    --min_identity 90 \
    --min_coverage 80

if [ $? -ne 0 ]; then
    echo "Error: Failed to extract coding sequences"
    exit 1
fi

# Count extracted sequences
if [ -d "output/operon_sequences" ]; then
    gene_count=$(ls output/operon_sequences/*.fasta 2>/dev/null | wc -l)
    echo "✓ Extracted sequences for $gene_count genes"
fi

# Also run BLAST-based analysis for comparison
echo "Running BLAST-based analysis for comparison..."
python extract_sequences_from_blast.py \
    --blast-csv ../03_blast_search/output/all_blast_hits_complete.csv \
    --output-dir output/diversity_analysis \
    --min-identity 90 \
    --min-coverage 80 \
    --analysis-only

echo "✓ Both sequence extraction and BLAST-based analysis completed"

################################################################################
# STEP 2: ANALYZE PROMOTER CONSERVATION FROM BLAST DATA
################################################################################
echo ""
echo "STEP 2: Analyzing promoter conservation from BLAST results..."
echo "-----------------------------------------------------------"
python create_promoter_msa_from_blast.py \
    --blast-dir ../03_blast_search/output/blast_results \
    --output-dir output/diversity_analysis

if [ $? -ne 0 ]; then
    echo "Warning: Promoter analysis failed, continuing..."
else
    echo "✓ Promoter conservation analysis completed"
fi

################################################################################
# STEP 3: CREATE MULTIPLE SEQUENCE ALIGNMENTS
################################################################################
echo ""
echo "STEP 3: Creating multiple sequence alignments..."
echo "----------------------------------------------"

# Create MSAs for coding sequences
echo "Creating MSAs for coding genes..."
python create_msa.py \
    --sequences-dir output/operon_sequences \
    --output-dir output/msa \
    --threads $THREADS

if [ $? -ne 0 ]; then
    echo "Error: Failed to create coding sequence alignments"
    exit 1
fi

# Count created alignments
coding_alignments=$(ls output/msa/dna_alignments/*_aligned.fasta 2>/dev/null | wc -l)
echo "✓ Created $coding_alignments coding sequence alignments"

################################################################################
# STEP 4: ANALYZE SEQUENCE DIVERSITY AND CREATE CONSERVATION PLOTS
################################################################################
echo ""
echo "STEP 4: Analyzing sequence diversity and creating conservation plots..."
echo "--------------------------------------------------------------------"

# Run diversity analysis with conservation plots
echo "Analyzing coding gene diversity..."
python analyze_diversity.py \
    --msa_dir output/msa \
    --output_dir output/diversity_analysis

if [ $? -ne 0 ]; then
    echo "Error: Failed to analyze coding gene diversity"
    exit 1
fi

echo "✓ Conservation plots and diversity analysis completed"

################################################################################
# STEP 5: CREATE GAP ANALYSIS PLOTS
################################################################################
echo ""
echo "STEP 5: Creating gap analysis plots..."
echo "------------------------------------"
python create_gap_analysis.py \
    --output-dir output/diversity_analysis

if [ $? -ne 0 ]; then
    echo "Warning: Gap analysis failed, continuing..."
else
    echo "✓ Gap analysis completed"
fi

################################################################################
# STEP 6: GENERATE FINAL SUMMARY
################################################################################
echo ""
echo "STEP 5: Generating final summary..."
echo "---------------------------------"

# Count final outputs
echo "Final Output Summary:"
echo "===================="

if [ -d "output/diversity_analysis" ]; then
    png_count=$(ls output/diversity_analysis/*.png 2>/dev/null | wc -l)
    csv_count=$(ls output/diversity_analysis/*.csv 2>/dev/null | wc -l)
    echo "  📊 Plots created: $png_count PNG files"
    echo "  📋 Data files: $csv_count CSV files"
    echo ""
    echo "Key output files:"
    echo "  - output/diversity_analysis/diversity_results.csv"
    echo "  - output/diversity_analysis/conservation_profiles.png"
    echo "  - output/diversity_analysis/promoter_conservation_profile.png"
    echo "  - output/diversity_analysis/gap_analysis_summary.png"
    echo "  - output/diversity_analysis/*_gap_profile.png"
fi

# Display diversity results summary
if [ -f "output/diversity_analysis/blast_based_diversity_results.csv" ]; then
    echo ""
    echo "BLAST-Based Diversity Results:"
    echo "============================="
    python -c "
import pandas as pd
import os

if os.path.exists('output/diversity_analysis/blast_based_diversity_results.csv'):
    df = pd.read_csv('output/diversity_analysis/blast_based_diversity_results.csv')
    print(f'Genes analyzed: {len(df)}')
    print(f'Average identity: {df[\"avg_identity\"].mean():.2f}%')
    print(f'Identity range: {df[\"min_identity\"].min():.2f}% - {df[\"max_identity\"].max():.2f}%')
    print(f'Average coverage: {df[\"avg_coverage\"].mean():.1f}%')
    
    print(f'\\nMost conserved gene: {df.loc[df[\"avg_identity\"].idxmax(), \"gene\"]} ({df[\"avg_identity\"].max():.2f}%)')
    print(f'Least conserved gene: {df.loc[df[\"avg_identity\"].idxmin(), \"gene\"]} ({df[\"avg_identity\"].min():.2f}%)')
"
fi

# Check promoter analysis
echo ""
echo "Promoter Analysis Summary:"
echo "========================"
python -c "
import pandas as pd
import os

blast_file = '../03_blast_search/output/all_blast_hits_complete.csv'
if os.path.exists(blast_file):
    df = pd.read_csv(blast_file)
    promoter_hits = df[df['element_name'] == 'promoter']
    if not promoter_hits.empty:
        best_hits = promoter_hits.loc[promoter_hits.groupby('genome_id')['pident'].idxmax()]
        print(f'Promoters analyzed: {len(best_hits)} genomes')
        print(f'Average promoter identity: {best_hits[\"pident\"].mean():.2f}%')
        print(f'Promoter identity range: {best_hits[\"pident\"].min():.1f}% - {best_hits[\"pident\"].max():.1f}%')
    else:
        print('No promoter data found')
else:
    print('BLAST results not found')
"

echo ""
echo "=========================================="
echo "COMPLETE DIVERSITY ANALYSIS FINISHED"
echo "=========================================="
echo "Completed at $(date)"
echo ""
echo "All results saved to: output/diversity_analysis/"
echo ""
echo "Next steps:"
echo "  1. View plots: output/diversity_analysis/*.png"
echo "  2. Check results: output/diversity_analysis/*.csv"
echo "  3. Read summary: output/diversity_analysis/diversity_summary.txt"