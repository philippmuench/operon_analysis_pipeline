#!/bin/bash
#SBATCH --job-name=core_gene_analysis
#SBATCH --output=core_analysis_%j.out
#SBATCH --error=core_analysis_%j.err
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=90
#SBATCH --mem=128G
#SBATCH --partition=cpu

# Parse command line arguments
START_STEP=1
HELP=false

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Core Gene Analysis Pipeline"
    echo ""
    echo "Options:"
    echo "  --start-step STEP    Start pipeline from specific step (1-5, default: 1)"
    echo "  --help               Show this help message"
    echo ""
    echo "Available steps:"
    echo "  1. Identify core genes (≥95% prevalence)"
    echo "  2. Extract sequences from Prokka output"
    echo "  3. Create MSAs for core genes"
    echo "  4. Calculate diversity metrics (optional)"
    echo "  5. Threshold analysis plots"
    echo ""
    echo "Examples:"
    echo "  $0                    # Run complete pipeline"
    echo "  $0 --start-step 3     # Start from MSA creation"
    echo "  $0 --start-step 5     # Only run threshold analysis"
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --start-step)
            START_STEP="$2"
            shift 2
            ;;
        --help)
            HELP=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

if [ "$HELP" = true ]; then
    usage
    exit 0
fi

# Validate start step
if ! [[ "$START_STEP" =~ ^[1-5]$ ]]; then
    echo "❌ Error: Invalid start step '$START_STEP'. Must be 1-5."
    usage
    exit 1
fi

# Set MAFFT temporary directory
export MAFFT_TMPDIR="/vol/tmp"

# Validate MAFFT_TMPDIR is properly set
echo "Validating MAFFT configuration..."
if [ -z "$MAFFT_TMPDIR" ]; then
    echo "❌ Error: MAFFT_TMPDIR environment variable is not set"
    echo "   This will cause MAFFT to use /tmp which may cause disk space issues"
    echo "   Please set: export MAFFT_TMPDIR=\"/vol/tmp\""
    exit 1
fi

if [ ! -d "$MAFFT_TMPDIR" ]; then
    echo "❌ Error: MAFFT_TMPDIR directory does not exist: $MAFFT_TMPDIR"
    echo "   Please create the directory or use a different path"
    exit 1
fi

if [ ! -w "$MAFFT_TMPDIR" ]; then
    echo "❌ Error: MAFFT_TMPDIR is not writable: $MAFFT_TMPDIR"
    echo "   Please check permissions on the temporary directory"
    exit 1
fi

echo "✅ MAFFT_TMPDIR validated: $MAFFT_TMPDIR"

# Core Gene Analysis Pipeline
# ===========================
# 1. Identify core genes (≥95% prevalence)
# 2. Extract sequences from Prokka output
# 3. Create MSAs for core genes
# 4. Calculate diversity metrics

echo "=================================================="
echo "Core Gene Analysis Pipeline"
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Starting from step: $START_STEP"
echo "=================================================="

# Change to script directory
cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/04_core_gene_analysis

# Step 1: Identify core genes
if [ $START_STEP -le 1 ]; then
    echo ""
    echo "Step 1: Identifying core genes..."
    echo "=================================="
    python identify_core_genes.py

    if [ $? -ne 0 ]; then
        echo "❌ Error: Core gene identification failed"
        exit 1
    fi

    # Check if core genes were found
    if [ ! -f "output/core_genes_95pct.txt" ]; then
        echo "❌ Error: Core genes file not created"
        exit 1
    fi

    CORE_GENE_COUNT=$(wc -l < output/core_genes_95pct.txt)
    echo "✅ Found $CORE_GENE_COUNT core genes"
else
    echo "⏭️  Skipping Step 1: Core gene identification"
    # Still need to count genes for reporting
    if [ -f "output/core_genes_95pct.txt" ]; then
        CORE_GENE_COUNT=$(wc -l < output/core_genes_95pct.txt)
        echo "   Found existing: $CORE_GENE_COUNT core genes"
    else
        echo "❌ Error: output/core_genes_95pct.txt not found. Run from step 1 first."
        exit 1
    fi
fi

# Step 2: Extract sequences
if [ $START_STEP -le 2 ]; then
    echo ""
    echo "Step 2: Extracting core gene sequences..."
    echo "=========================================="
    python extract_core_sequences.py --threads $SLURM_CPUS_PER_TASK

    if [ $? -ne 0 ]; then
        echo "❌ Error: Sequence extraction failed"
        exit 1
    fi

    # Check sequences were extracted
    if [ ! -d "output/core_gene_sequences" ]; then
        echo "❌ Error: Sequence directory not created"
        exit 1
    fi

    FASTA_COUNT=$(find output/core_gene_sequences -name "*.fasta" | wc -l)
    echo "✅ Extracted sequences for $FASTA_COUNT genes"
else
    echo "⏭️  Skipping Step 2: Sequence extraction"
    # Still need to count sequences for reporting
    if [ -d "output/core_gene_sequences" ]; then
        FASTA_COUNT=$(find output/core_gene_sequences -name "*.fasta" | wc -l)
        echo "   Found existing: $FASTA_COUNT sequence files"
    else
        echo "❌ Error: output/core_gene_sequences/ not found. Run from step 2 first."
        exit 1
    fi
fi

# Step 3: Create MSAs
if [ $START_STEP -le 3 ]; then
    echo ""
    echo "Step 3: Creating multiple sequence alignments..."
    echo "================================================="
    python create_core_gene_msa.py --threads $SLURM_CPUS_PER_TASK

    if [ $? -ne 0 ]; then
        echo "❌ Error: MSA creation failed"
        exit 1
    fi

    # Check alignments were created
    if [ ! -d "output/core_gene_alignments" ]; then
        echo "❌ Error: Alignment directory not created"
        exit 1
    fi

    ALIGNMENT_COUNT=$(find output/core_gene_alignments -name "*_aligned.fasta" | wc -l)
    echo "✅ Created $ALIGNMENT_COUNT alignments"
else
    echo "⏭️  Skipping Step 3: MSA creation"
    # Still need to count alignments for reporting
    if [ -d "output/core_gene_alignments" ]; then
        ALIGNMENT_COUNT=$(find output/core_gene_alignments -name "*_aligned.fasta" | wc -l)
        echo "   Found existing: $ALIGNMENT_COUNT alignment files"
    else
        echo "❌ Error: output/core_gene_alignments/ not found. Run from step 3 first."
        exit 1
    fi
fi

# Step 4: Calculate diversity metrics (optional - create if needed)
if [ $START_STEP -le 4 ] && [ -f "calculate_diversity.py" ]; then
    echo ""
    echo "Step 4: Calculating diversity metrics..."
    echo "========================================"
    python calculate_diversity.py --threads $SLURM_CPUS_PER_TASK --progress-every 20 --checkpoint-every 200
    
    if [ $? -eq 0 ]; then
        echo "✅ Diversity analysis complete"
    else
        echo "⚠️  Warning: Diversity analysis failed"
    fi
elif [ $START_STEP -le 4 ]; then
    echo "⏭️  Skipping Step 4: Diversity metrics (calculate_diversity.py not found)"
else
    echo "⏭️  Skipping Step 4: Diversity metrics"
fi

# Summary
echo ""
echo "=================================================="
echo "Core Gene Analysis Complete!"
echo "Finished: $(date)"
echo "=================================================="
echo "📊 Results:"
echo "   - Core genes identified: $CORE_GENE_COUNT"
echo "   - Sequences extracted: $FASTA_COUNT genes"
echo "   - MSAs created: $ALIGNMENT_COUNT"
echo ""
# Step 5: Threshold Analysis
if [ $START_STEP -le 5 ]; then
    echo ""
    echo "Step 5: Creating threshold analysis plots..."
    echo "=========================================="
    python plot_core_gene_thresholds.py \
        --core-genes-file output/gene_prevalence_stats.csv \
        --prokka-dir ../01_prokka_annotation/output/prokka_results \
        --output-dir output

    if [ $? -eq 0 ]; then
        echo "✅ Threshold analysis complete"
    else
        echo "⚠️  Warning: Threshold analysis failed"
    fi
else
    echo "⏭️  Skipping Step 5: Threshold analysis"
fi


echo "📁 Output directories:"
echo "   - output/core_genes_95pct.txt"
echo "   - output/gene_prevalence_stats.csv"
echo "   - output/core_gene_sequences/"
echo "   - output/core_gene_alignments/"
echo "   - output/core_gene_threshold_*.png"
echo ""

# Check final results
if [ $ALIGNMENT_COUNT -gt 0 ]; then
    echo "✅ Pipeline completed successfully!"
    echo ""
    echo "🧬 Ready for diversity analysis:"
    echo "   - Use MSAs in output/core_gene_alignments/"
    echo "   - Compare with operon gene diversity"
    echo "   - Calculate conservation metrics"
    exit 0
else
    echo "❌ Pipeline failed - no alignments created"
    exit 1
fi