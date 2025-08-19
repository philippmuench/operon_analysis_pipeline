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
# STEP 1: VERIFY EXISTING MSA DATA FROM PREVIOUS STEPS
################################################################################
echo ""
echo "STEP 1: Verifying existing MSA data from previous steps..."
echo "--------------------------------------------------------"

# Check Strategy A (Prokka-based) MSAs
STRATEGY_A_DIR="../05_operon_assembly_extraction/output/mappings/aa_nt_mapping/prokka/msa/dna_alignments"
STRATEGY_D_DIR="../05_operon_assembly_extraction/output/mappings/aa_nt_mapping/assemblies/msa/dna_alignments"
CORE_GENE_DIR="../04_core_gene_analysis/output/core_gene_alignments"

echo "Checking Strategy A (Prokka-based) MSAs..."
if [ -d "$STRATEGY_A_DIR" ]; then
    strategy_a_count=$(ls "$STRATEGY_A_DIR"/*_aligned.fasta 2>/dev/null | wc -l)
    echo "✓ Found $strategy_a_count Strategy A operon gene MSAs"
else
    echo "❌ Strategy A MSAs not found at $STRATEGY_A_DIR"
    exit 1
fi

echo "Checking Strategy D (Assembly-based) MSAs..."
if [ -d "$STRATEGY_D_DIR" ]; then
    strategy_d_count=$(ls "$STRATEGY_D_DIR"/*_aligned.fasta 2>/dev/null | wc -l)
    echo "✓ Found $strategy_d_count Strategy D operon gene MSAs"
else
    echo "⚠️  Strategy D MSAs not found at $STRATEGY_D_DIR - may need to complete Step 05"
    echo "   Continuing with Strategy A only..."
    STRATEGY_D_DIR=""
fi

echo "Checking core gene MSAs..."
if [ -d "$CORE_GENE_DIR" ]; then
    core_gene_count=$(ls "$CORE_GENE_DIR"/*.fasta 2>/dev/null | wc -l)
    echo "✓ Found $core_gene_count core gene MSAs"
else
    echo "❌ Core gene MSAs not found at $CORE_GENE_DIR"
    exit 1
fi

# Create symbolic links to existing data for compatibility
echo ""
echo "Creating symbolic links to existing MSA data..."
mkdir -p output/msa_strategy_a/dna_alignments
mkdir -p output/msa_strategy_d/dna_alignments
mkdir -p output/core_gene_msa

# Link Strategy A MSAs
ln -sf "$(realpath $STRATEGY_A_DIR)"/* output/msa_strategy_a/dna_alignments/ 2>/dev/null || true

# Link Strategy D MSAs if available
if [ -n "$STRATEGY_D_DIR" ] && [ -d "$STRATEGY_D_DIR" ]; then
    ln -sf "$(realpath $STRATEGY_D_DIR)"/* output/msa_strategy_d/dna_alignments/ 2>/dev/null || true
fi

# Link core gene MSAs
ln -sf "$(realpath $CORE_GENE_DIR)"/* output/core_gene_msa/ 2>/dev/null || true

echo "✓ Data preparation completed using existing high-quality MSAs"

################################################################################
# STEP 2: ANALYZE PROMOTER CONSERVATION FROM EXISTING DATA
################################################################################
echo ""
echo "STEP 2: Using promoter conservation from existing operon extraction..."
echo "--------------------------------------------------------------------"

# Use existing promoter MSAs from Step 05
PROMOTER_A_DIR="../05_operon_assembly_extraction/output/mappings/aa_nt_mapping/prokka/msa/noncoding_alignments"
PROMOTER_D_DIR="../05_operon_assembly_extraction/output/mappings/aa_nt_mapping/assemblies/msa/noncoding_alignments"

mkdir -p output/promoter_msa

if [ -d "$PROMOTER_A_DIR" ]; then
    ln -sf "$(realpath $PROMOTER_A_DIR)"/* output/promoter_msa/ 2>/dev/null || true
    echo "✓ Linked Strategy A promoter MSAs"
elif [ -d "$PROMOTER_D_DIR" ]; then
    ln -sf "$(realpath $PROMOTER_D_DIR)"/* output/promoter_msa/ 2>/dev/null || true
    echo "✓ Linked Strategy D promoter MSAs"
else
    echo "⚠️  No promoter MSAs found - running BLAST-based analysis as fallback..."
    python create_promoter_msa_from_blast.py \
        --blast-dir ../03_blast_search/output/blast_results \
        --output-dir output/diversity_analysis
fi

################################################################################
# STEP 3: COMPREHENSIVE DIVERSITY ANALYSIS (STRATEGIES A, D, AND CORE GENES)
################################################################################
echo ""
echo "STEP 3: Comprehensive diversity analysis..."
echo "=========================================="

# Strategy A: Prokka-based operon gene analysis
echo ""
echo "Analyzing Strategy A (Prokka-based) operon genes..."
echo "--------------------------------------------------"
python analyze_diversity.py \
    --msa_dir output/msa_strategy_a \
    --output_dir output/diversity_analysis_strategy_a

if [ $? -eq 0 ]; then
    echo "✓ Strategy A diversity analysis completed"
else
    echo "⚠️  Strategy A analysis failed, continuing..."
fi

# Strategy D: Assembly-based operon gene analysis (if available)
if [ -n "$STRATEGY_D_DIR" ] && [ -d "$STRATEGY_D_DIR" ]; then
    echo ""
    echo "Analyzing Strategy D (Assembly-based) operon genes..."
    echo "---------------------------------------------------"
    python analyze_diversity.py \
        --msa_dir output/msa_strategy_d \
        --output_dir output/diversity_analysis_strategy_d

    if [ $? -eq 0 ]; then
        echo "✓ Strategy D diversity analysis completed"
    else
        echo "⚠️  Strategy D analysis failed, continuing..."
    fi
else
    echo ""
    echo "⏭️  Skipping Strategy D analysis (data not available)"
fi

# Core genes analysis
echo ""
echo "Analyzing core genes for baseline comparison..."
echo "---------------------------------------------"
python analyze_all_core_genes.py \
    --core_dir output/core_gene_msa \
    --output_dir output/core_gene_analysis

if [ $? -eq 0 ]; then
    echo "✓ Core gene diversity analysis completed"
else
    echo "⚠️  Core gene analysis failed, continuing..."
fi

echo ""
echo "✓ Comprehensive diversity analysis completed"

################################################################################
# STEP 4: CREATE COMPARATIVE ANALYSIS AND FINAL SUMMARY
################################################################################
echo ""
echo "STEP 4: Creating comparative analysis and final summary..."
echo "--------------------------------------------------------"

# Generate comparative diversity summary
echo "Creating comparative diversity summary..."
python generate_diversity_summary.py \
    --strategy_a_dir output/diversity_analysis_strategy_a \
    --strategy_d_dir output/diversity_analysis_strategy_d \
    --core_gene_dir output/core_gene_analysis \
    --output_dir output/comparative_analysis

if [ $? -eq 0 ]; then
    echo "✓ Comparative analysis completed"
else
    echo "⚠️  Comparative analysis failed, showing individual results..."
fi

# Create gap analysis for all strategies
echo ""
echo "Creating gap analysis for all datasets..."
python create_gap_analysis.py \
    --output-dir output/comparative_analysis

################################################################################
# STEP 5: GENERATE COMPREHENSIVE FINAL SUMMARY
################################################################################
echo ""
echo "STEP 5: Generating comprehensive final summary..."
echo "-----------------------------------------------"

# Count final outputs
echo ""
echo "COMPREHENSIVE OUTPUT SUMMARY:"
echo "============================"

# Strategy A outputs
if [ -d "output/diversity_analysis_strategy_a" ]; then
    strategy_a_plots=$(ls output/diversity_analysis_strategy_a/*.png 2>/dev/null | wc -l)
    strategy_a_data=$(ls output/diversity_analysis_strategy_a/*.csv 2>/dev/null | wc -l)
    echo "📊 Strategy A (Prokka-based): $strategy_a_plots plots, $strategy_a_data data files"
fi

# Strategy D outputs
if [ -d "output/diversity_analysis_strategy_d" ]; then
    strategy_d_plots=$(ls output/diversity_analysis_strategy_d/*.png 2>/dev/null | wc -l)
    strategy_d_data=$(ls output/diversity_analysis_strategy_d/*.csv 2>/dev/null | wc -l)
    echo "📊 Strategy D (Assembly-based): $strategy_d_plots plots, $strategy_d_data data files"
fi

# Core gene outputs
if [ -d "output/core_gene_analysis" ]; then
    core_plots=$(ls output/core_gene_analysis/*.png 2>/dev/null | wc -l)
    core_data=$(ls output/core_gene_analysis/*.csv 2>/dev/null | wc -l)
    echo "📊 Core genes baseline: $core_plots plots, $core_data data files"
fi

# Comparative outputs
if [ -d "output/comparative_analysis" ]; then
    comp_plots=$(ls output/comparative_analysis/*.png 2>/dev/null | wc -l)
    comp_data=$(ls output/comparative_analysis/*.csv 2>/dev/null | wc -l)
    echo "📊 Comparative analysis: $comp_plots plots, $comp_data data files"
fi

echo ""
echo "KEY OUTPUT DIRECTORIES:"
echo "======================"
echo "  📁 output/diversity_analysis_strategy_a/  (Prokka-based operon analysis)"
echo "  📁 output/diversity_analysis_strategy_d/  (Assembly-based operon analysis)"
echo "  📁 output/core_gene_analysis/             (Core genes baseline)"
echo "  📁 output/comparative_analysis/           (Strategy comparison)"
echo "  📁 output/promoter_msa/                   (Promoter sequences)"

echo ""
echo "=========================================="
echo "COMPLETE DIVERSITY ANALYSIS FINISHED"
echo "=========================================="
echo "Completed at $(date)"
echo ""
echo "NEXT STEPS FOR ANALYSIS:"
echo "======================="
echo "  1. Compare Strategy A vs D results for operon genes"
echo "  2. Compare operon conservation vs core gene baseline"
echo "  3. Examine promoter conservation patterns"
echo "  4. Proceed to dN/dS analysis using selected strategy"
echo ""
echo "KEY RESULT FILES:"
echo "================"
echo "  📊 output/diversity_analysis_strategy_a/diversity_results.csv"
echo "  📊 output/diversity_analysis_strategy_d/diversity_results.csv"
echo "  📊 output/core_gene_analysis/core_gene_diversity_results.csv"
echo "  📊 output/comparative_analysis/strategy_comparison.csv"
echo ""
echo "VISUALIZATION FILES:"
echo "==================="
echo "  🎯 output/*/conservation_profiles.png"
echo "  🎯 output/comparative_analysis/diversity_comparison.png"
echo "  🎯 output/promoter_msa/promoter_conservation.png"