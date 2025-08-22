#!/bin/bash
#SBATCH --job-name=operon_extraction
#SBATCH --output=operon_extraction_%j.out
#SBATCH --error=operon_extraction_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=cpu
#SBATCH --qos=verylong

# Operon Sequence Extraction and MSA Pipeline - Clean Version
# ============================================================
# Primary analysis: Strategy D (assembly-based tblastn)
# Optional: Gene boundary analysis (Strategy A - prokka-based)

# Parse command line arguments
START_STEP=1
HELP=false
WITH_GENE_BOUNDARY=false
METRICS_ONLY=false

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Operon Sequence Extraction Pipeline (Clean Version)"
    echo ""
    echo "Options:"
    echo "  --start-step STEP        Start pipeline from specific step (1-5, default: 1)"
    echo "  --with-gene-boundary     Include gene boundary analysis (Strategy A)"
    echo "  --metrics-only           Only calculate conservation metrics (Step 4)"
    echo "  --help                   Show this help message"
    echo ""
    echo "Available steps:"
    echo "  1. Extract operon gene sequences from assemblies (Strategy D)"
    echo "  2. Create MSAs from gene sequences (DNA & protein)"
    echo "  3. Extract promoter sequences and create promoter MSA"
    echo "  4. Calculate conservation metrics (plots in step 06)"
    echo "  5. Generate comprehensive summary and statistics"
    echo ""
    echo "Examples:"
    echo "  $0                              # Run complete pipeline"
    echo "  $0 --with-gene-boundary         # Include gene boundary analysis"
    echo "  $0 --start-step 4               # Resume from step 4"
    echo "  $0 --metrics-only               # Quick re-run of metrics calculation"
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --start-step)
            START_STEP="$2"
            shift 2
            ;;
        --with-gene-boundary)
            WITH_GENE_BOUNDARY=true
            shift
            ;;
        --metrics-only)
            METRICS_ONLY=true
            START_STEP=4
            shift
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

echo "=================================================="
if [ "$METRICS_ONLY" = true ]; then
    echo "Operon Conservation Metrics Calculation"
    echo "(Metrics-only mode - using existing MSAs)"
else
    echo "Operon Sequence Extraction Pipeline (Clean)"
fi
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Starting from step: $START_STEP"
if [ "$WITH_GENE_BOUNDARY" = true ]; then
    echo "Including: Gene boundary analysis"
fi
echo "=================================================="

# Change to script directory
cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/05_operon_assembly_extraction

# Ensure MAFFT temporary directory uses fast local scratch by default
MAFFT_TMP_DEFAULT="/vol/tmp"
mkdir -p "${MAFFT_TMPDIR:-$MAFFT_TMP_DEFAULT}" || true
export MAFFT_TMPDIR="${MAFFT_TMPDIR:-$MAFFT_TMP_DEFAULT}"
export TMPDIR="${TMPDIR:-$MAFFT_TMPDIR}"
echo "Using MAFFT_TMPDIR=$MAFFT_TMPDIR (TMPDIR=$TMPDIR)"

# Create main output directories
echo ""
echo "Setting up output directories..."
echo "================================"
mkdir -p output/{sequences,msa,plots,noncoding_sequences}
mkdir -p output/msa/{dna_alignments,noncoding_alignments}

# Setup gene boundary analysis directory if requested
if [ "$WITH_GENE_BOUNDARY" = true ]; then
    mkdir -p output/gene_boundary_analysis/{sequences,msa,plots}
    mkdir -p output/gene_boundary_analysis/msa/{dna_alignments,noncoding_alignments}
fi

# Initialize conda environment
echo ""
echo "Activating conda environment..."
eval "$(/home/pmuench/miniconda3/bin/conda shell.bash hook)"
conda activate efs_diversity
echo "Conda environment activated: $CONDA_DEFAULT_ENV"

THREADS=${SLURM_CPUS_PER_TASK:-8}

########## Step 1: Extract sequences from assemblies (Strategy D) ##########
if [ $START_STEP -le 1 ] && [ "$METRICS_ONLY" = false ]; then
    echo ""
    echo "Step 1: Extracting operon gene sequences from assemblies..."
    echo "==========================================================="
    echo "🧬 Using Strategy D: tblastn (aa→nt) with raw assemblies"
    echo "   This is the primary analysis for phylogenetic studies"
    echo ""
    
    python operon_pipeline.py extract-sequences \
        --blast-dir ../03_blast_search/output/blast_results \
        --output-dir output/sequences \
        --assemblies-dir ../../Efs_assemblies \
        --min-identity 90 --min-coverage 80 --source assemblies
    
    if [ $? -eq 0 ]; then
        echo "✅ Sequence extraction completed"
        
        # Verify extraction produced output
        test -s "output/sequences/frpC.fasta" || { echo "❌ Error: No sequences extracted"; exit 2; }
        
        # Count extracted sequences
        TOTAL_FILES=$(find output/sequences -name "*.fasta" | wc -l)
        for gene in frpC glpC ptsD ptsC ptsB ptsA fruR; do
            if [ -f "output/sequences/${gene}.fasta" ]; then
                COUNT=$(grep -c "^>" "output/sequences/${gene}.fasta")
                echo "   ${gene}: $COUNT sequences"
            fi
        done
    else
        echo "❌ Error: Sequence extraction failed"
        exit 1
    fi
    
    # Optional: Gene boundary analysis
    if [ "$WITH_GENE_BOUNDARY" = true ]; then
        echo ""
        echo "Extracting sequences for gene boundary analysis..."
        echo "   Using Strategy A: tblastn with Prokka annotations"
        
        python operon_pipeline.py extract-sequences \
            --blast-dir ../03_blast_search/output/blast_results \
            --output-dir output/gene_boundary_analysis/sequences \
            --prokka-dir ../01_prokka_annotation/output/prokka_results \
            --min-identity 90 --min-coverage 80 --source prokka
        
        if [ $? -eq 0 ]; then
            echo "✅ Gene boundary extraction completed"
        else
            echo "⚠️  Warning: Gene boundary extraction failed"
        fi
    fi
else
    echo "⏭️  Skipping Step 1: Sequence extraction"
fi

########## Step 2: Create MSAs ##########
if [ $START_STEP -le 2 ] && [ "$METRICS_ONLY" = false ]; then
    echo ""
    echo "Step 2: Creating multiple sequence alignments..."
    echo "================================================"
    echo "🧬 Aligning extracted sequences with MAFFT"
    echo ""
    
    python operon_pipeline.py create-msa \
        --sequences-dir output/sequences \
        --output-dir output/msa \
        --threads "$THREADS"
    
    if [ $? -eq 0 ]; then
        echo "✅ MSA creation completed"
        
        # Count alignments
        DNA_COUNT=$(find output/msa/dna_alignments -name "*_aligned.fasta" 2>/dev/null | wc -l)
        echo "   Created $DNA_COUNT DNA alignments"
    else
        echo "❌ Error: MSA creation failed"
        exit 1
    fi
    
    # Optional: Gene boundary MSAs
    if [ "$WITH_GENE_BOUNDARY" = true ] && [ -d "output/gene_boundary_analysis/sequences" ]; then
        echo ""
        echo "Creating MSAs for gene boundary analysis..."
        
        python operon_pipeline.py create-msa \
            --sequences-dir output/gene_boundary_analysis/sequences \
            --output-dir output/gene_boundary_analysis/msa \
            --threads "$THREADS"
        
        if [ $? -eq 0 ]; then
            echo "✅ Gene boundary MSA creation completed"
        else
            echo "⚠️  Warning: Gene boundary MSA creation failed"
        fi
    fi
else
    echo "⏭️  Skipping Step 2: MSA creation"
fi

########## Step 3: Extract promoter sequences ##########
if [ $START_STEP -le 3 ] && [ "$METRICS_ONLY" = false ]; then
    echo ""
    echo "Step 3: Extracting and aligning promoter sequences..."
    echo "====================================================="
    echo "🧬 Processing non-coding regulatory regions"
    echo ""
    
    python operon_pipeline.py extract-noncoding \
        --blast-dir ../03_blast_search/output/blast_results \
        --assemblies-dir ../../Efs_assemblies \
        --output-dir output/noncoding_sequences \
        --min-identity 80 --min-coverage 70 --source assemblies
    
    if [ $? -eq 0 ]; then
        echo "✅ Promoter extraction completed"
        
        # Create promoter MSA
        if [ -f "output/noncoding_sequences/promoter.fasta" ]; then
            echo "Creating promoter MSA..."
            mafft --auto --thread "$THREADS" output/noncoding_sequences/promoter.fasta \
                > output/msa/noncoding_alignments/promoter_aligned.fasta 2>/dev/null
            
            if [ $? -eq 0 ]; then
                PROMOTER_COUNT=$(grep -c "^>" output/msa/noncoding_alignments/promoter_aligned.fasta)
                echo "✅ Promoter MSA created with $PROMOTER_COUNT sequences"
            fi
        fi
    else
        echo "⚠️  Warning: Promoter extraction failed"
    fi
else
    echo "⏭️  Skipping Step 3: Promoter extraction"
fi

########## Step 4: Create conservation metrics (plots moved to step 06) ##########
if [ $START_STEP -le 4 ]; then
    echo ""
    echo "Step 4: Calculating conservation metrics..."
    echo "=================================================="
    if [ "$METRICS_ONLY" = true ]; then
        echo "📊 Running in metrics-only mode (using existing MSAs)"
    else
        echo "📊 Generating conservation scores (plots will be created in step 06)"
    fi
    echo ""
    
    # Generate conservation metrics CSV only (no plotting)
    if [ -d "output/msa/dna_alignments" ]; then
        # Count alignment files first
        ALIGN_COUNT=$(find output/msa/dna_alignments -name "*_aligned.fasta" 2>/dev/null | wc -l)
        echo "Found $ALIGN_COUNT alignment files in output/msa/dna_alignments"
        
        if [ $ALIGN_COUNT -eq 0 ]; then
            echo "❌ Error: No alignment files found"
            echo "   Please ensure Step 2 (MSA creation) completed successfully"
            exit 1
        fi
        
        echo "Calculating conservation metrics from alignments..."
        
        # Generate basic metrics (for backward compatibility)
        python operon_pipeline.py create-summary \
            --msa-dir output/msa \
            --output-file output/operon_conservation_metrics.csv \
            --strategy-name "Primary Analysis"
        
        if [ $? -eq 0 ]; then
            echo "✅ Basic conservation metrics saved to output/operon_conservation_metrics.csv"
            
            # Show summary statistics
            echo ""
            echo "Conservation summary:"
            echo "--------------------"
            head -n 1 output/operon_conservation_metrics.csv
            tail -n +2 output/operon_conservation_metrics.csv | head -n 10
            echo "..."
            TOTAL_GENES=$(tail -n +2 output/operon_conservation_metrics.csv | wc -l)
            echo "Total genes analyzed: $TOTAL_GENES"
            
            # Calculate average conservation
            if [ $TOTAL_GENES -gt 0 ]; then
                AVG_CONS=$(tail -n +2 output/operon_conservation_metrics.csv | \
                           awk -F',' '{sum+=$4; count++} END {if(count>0) printf "%.3f", sum/count}')
                echo "Average conservation score: $AVG_CONS"
            fi
        else
            echo "❌ Error: Failed to generate basic conservation metrics"
            exit 1
        fi
        
        # Generate detailed metrics matching core gene format
        echo ""
        echo "Calculating detailed conservation metrics (matching core gene format)..."
        python operon_pipeline.py create-detailed-summary \
            --msa-dir output/msa \
            --output-file output/operon_conservation_metrics_detailed.csv
        
        if [ $? -eq 0 ]; then
            echo "✅ Detailed conservation metrics saved to output/operon_conservation_metrics_detailed.csv"
            echo ""
            echo "Detailed metrics preview:"
            echo "------------------------"
            head -n 5 output/operon_conservation_metrics_detailed.csv
        else
            echo "⚠️  Warning: Failed to generate detailed conservation metrics"
        fi
        
        # Generate gap fraction plots
        echo ""
        echo "Creating gap fraction plots..."
        python operon_pipeline.py create-gap-plots \
            --msa-dir output/msa/dna_alignments \
            --output-dir output/plots/gaps
        
        if [ $? -eq 0 ]; then
            echo "✅ Gap fraction plots saved to output/plots/gaps/"
            GAP_PLOT_COUNT=$(ls output/plots/gaps/*_gap_fraction_analysis.png 2>/dev/null | wc -l)
            echo "   Generated $GAP_PLOT_COUNT gap analysis plots"
        else
            echo "⚠️  Warning: Failed to generate gap fraction plots"
        fi
    else
        echo "❌ Error: No MSA alignments found in output/msa/dna_alignments"
        echo "   Please ensure Step 2 (MSA creation) completed successfully"
        exit 1
    fi
    
    # Optional: Gene boundary metrics
    if [ "$WITH_GENE_BOUNDARY" = true ] && [ -d "output/gene_boundary_analysis/msa/dna_alignments" ]; then
        echo ""
        echo "Calculating gene boundary analysis metrics..."
        
        python operon_pipeline.py create-summary \
            --msa-dir output/gene_boundary_analysis/msa \
            --output-file output/gene_boundary_analysis/gene_boundary_conservation_metrics.csv \
            --strategy-name "Gene Boundary Analysis"
        
        if [ $? -eq 0 ]; then
            echo "✅ Gene boundary metrics saved"
        fi
    fi
    
    echo ""
    echo "NOTE: To generate conservation plots, run:"
    echo "  cd ../06_diversity_analysis"
    echo "  sbatch run_analysis.sh --plots"
else
    echo "⏭️  Skipping Step 4: Conservation metrics"
fi

########## Step 5: Generate comprehensive summary ##########
if [ $START_STEP -le 5 ] && [ "$METRICS_ONLY" = false ]; then
    echo ""
    echo "Step 5: Generating comprehensive summary..."
    echo "==========================================="
    echo "📊 Creating extraction pipeline summary"
    echo ""
    
    # Create summary file
    SUMMARY_FILE="output/extraction_pipeline_summary.txt"
    
    echo "================================================" > "$SUMMARY_FILE"
    echo "Operon Extraction Pipeline Summary" >> "$SUMMARY_FILE"
    echo "Generated: $(date)" >> "$SUMMARY_FILE"
    echo "================================================" >> "$SUMMARY_FILE"
    echo "" >> "$SUMMARY_FILE"
    
    echo "PRIMARY ANALYSIS (Assembly-based extraction)" >> "$SUMMARY_FILE"
    echo "============================================" >> "$SUMMARY_FILE"
    
    # Count sequences
    for gene in frpC glpC ptsD ptsC ptsB ptsA fruR; do
        if [ -f "output/sequences/${gene}.fasta" ]; then
            COUNT=$(grep -c "^>" "output/sequences/${gene}.fasta")
            echo "  ${gene}: $COUNT sequences extracted" >> "$SUMMARY_FILE"
        fi
    done
    
    echo "" >> "$SUMMARY_FILE"
    echo "ALIGNMENTS" >> "$SUMMARY_FILE"
    echo "==========" >> "$SUMMARY_FILE"
    
    # Count alignments
    if [ -d "output/msa/dna_alignments" ]; then
        for align in output/msa/dna_alignments/*_aligned.fasta; do
            if [ -f "$align" ]; then
                BASENAME=$(basename "$align" _aligned.fasta)
                COUNT=$(grep -c "^>" "$align")
                LENGTH=$(awk '/^>/{if (seq){print length(seq); exit} next}{seq=seq""$0} END{print length(seq)}' "$align")
                echo "  ${BASENAME}: $COUNT sequences, ${LENGTH} bp alignment" >> "$SUMMARY_FILE"
            fi
        done
    fi
    
    if [ -f "output/msa/noncoding_alignments/promoter_aligned.fasta" ]; then
        PROM_COUNT=$(grep -c "^>" output/msa/noncoding_alignments/promoter_aligned.fasta)
        echo "  Promoter: $PROM_COUNT sequences aligned" >> "$SUMMARY_FILE"
    fi
    
    echo "" >> "$SUMMARY_FILE"
    echo "CONSERVATION METRICS" >> "$SUMMARY_FILE"
    echo "===================" >> "$SUMMARY_FILE"
    
    if [ -f "output/operon_conservation_metrics.csv" ]; then
        echo "  Metrics saved in: operon_conservation_metrics.csv" >> "$SUMMARY_FILE"
        
        # Show top conservation scores
        echo "" >> "$SUMMARY_FILE"
        echo "  Top conservation scores:" >> "$SUMMARY_FILE"
        tail -n +2 output/operon_conservation_metrics.csv | \
            sort -t',' -k4 -rn | head -n 3 | \
            while IFS=',' read -r gene seqs pos score gaps pairwise; do
                echo "    ${gene}: ${score} (conservation score)" >> "$SUMMARY_FILE"
            done
    fi
    
    echo "" >> "$SUMMARY_FILE"
    echo "OUTPUT FILES" >> "$SUMMARY_FILE"
    echo "============" >> "$SUMMARY_FILE"
    echo "  Sequences: output/sequences/*.fasta" >> "$SUMMARY_FILE"
    echo "  DNA alignments: output/msa/dna_alignments/*_aligned.fasta" >> "$SUMMARY_FILE"
    echo "  Conservation plots: output/plots/*.png" >> "$SUMMARY_FILE"
    echo "  Metrics: output/operon_conservation_metrics.csv" >> "$SUMMARY_FILE"
    
    if [ "$WITH_GENE_BOUNDARY" = true ]; then
        echo "" >> "$SUMMARY_FILE"
        echo "GENE BOUNDARY ANALYSIS" >> "$SUMMARY_FILE"
        echo "=====================" >> "$SUMMARY_FILE"
        echo "  Directory: output/gene_boundary_analysis/" >> "$SUMMARY_FILE"
        echo "  Purpose: Validation of gene boundaries using Prokka annotations" >> "$SUMMARY_FILE"
    fi
    
    echo "" >> "$SUMMARY_FILE"
    echo "================================================" >> "$SUMMARY_FILE"
    echo "Pipeline completed successfully at $(date)" >> "$SUMMARY_FILE"
    echo "================================================" >> "$SUMMARY_FILE"
    
    echo "✅ Summary saved to $SUMMARY_FILE"
    
    # Display summary
    cat "$SUMMARY_FILE"
else
    echo "⏭️  Skipping Step 5: Summary generation"
fi

echo ""
echo "=================================================="
echo "Pipeline completed at $(date)"
echo "=================================================="
echo ""
echo "📁 Main output directory: output/"
echo "   - Sequences: output/sequences/"
echo "   - Alignments: output/msa/"
echo "   - Metrics: output/operon_conservation_metrics.csv"
echo ""
echo "📊 To generate conservation plots:"
echo "   cd ../06_diversity_analysis && sbatch run_analysis.sh --plots"

if [ "$WITH_GENE_BOUNDARY" = true ]; then
    echo ""
    echo "📁 Gene boundary analysis: output/gene_boundary_analysis/"
fi

echo ""
echo "✨ Analysis complete!"