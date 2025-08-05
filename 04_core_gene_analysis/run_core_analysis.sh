#!/bin/bash
#SBATCH --job-name=core_gene_analysis
#SBATCH --output=core_analysis_%j.out
#SBATCH --error=core_analysis_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=30
#SBATCH --mem=64G
#SBATCH --partition=cpu

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
echo "=================================================="

# Change to script directory
cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/04_core_gene_analysis

# Step 1: Identify core genes
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

# Step 2: Extract sequences
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

# Step 3: Create MSAs
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

# Step 4: Calculate diversity metrics (optional - create if needed)
if [ -f "calculate_diversity.py" ]; then
    echo ""
    echo "Step 4: Calculating diversity metrics..."
    echo "========================================"
    python calculate_diversity.py
    
    if [ $? -eq 0 ]; then
        echo "✅ Diversity analysis complete"
    else
        echo "⚠️  Warning: Diversity analysis failed"
    fi
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