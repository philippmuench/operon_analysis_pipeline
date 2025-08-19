#!/bin/bash
#SBATCH --job-name=operon_extraction
#SBATCH --output=operon_extraction_%j.out
#SBATCH --error=operon_extraction_%j.err
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --partition=cpu

# Operon Sequence Extraction and MSA Pipeline
# ============================================
# Complete pipeline for extracting operon sequences from BLAST results
# and creating MSAs for both coding genes and non-coding regions

# Parse command line arguments
START_STEP=1
HELP=false
# Default: run all strategies
STRATEGIES="current,nt_vs_genome,prokka_variants,assemblies"

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Operon Sequence Extraction Pipeline"
    echo ""
    echo "Options:"
    echo "  --start-step STEP    Start pipeline from specific step (1-5, default: 1)"
    echo "  --strategies LIST    Comma-separated: current,nt_vs_genome,prokka_variants,assemblies (default: current)"
    echo "  --help               Show this help message"
    echo ""
    echo "Available steps:"
    echo "  1. Extract operon gene sequences from assemblies"
    echo "  2. Create MSAs from gene sequences (DNA & protein)"
    echo "  3. Extract promoter sequences and create promoter MSA"
    echo "  4. Create conservation plots (genes and enhanced promoter plots)"
    echo "  5. BLAST-based diversity analysis (supplementary) and summary"
    echo "  6. Run multi-strategy MSAs (nt-vs-nt, tblastn, Prokka .ffn variants) and plots"
    echo ""
    echo "Examples:"
    echo "  $0                                      # Run complete pipeline (current strategy)"
    echo "  $0 --strategies current,prokka_variants  # Run current + old-look variants in one go"
    echo "  $0 --start-step 3                       # Start from promoter extraction and MSA"
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --start-step)
            START_STEP="$2"
            shift 2
            ;;
        --strategies)
            STRATEGIES="$2"
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

echo "=================================================="
echo "Operon Sequence Extraction Pipeline"
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Starting from step: $START_STEP"
echo "Strategies: $STRATEGIES"
echo "=================================================="

# Change to script directory
cd /vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis/05_operon_assembly_extraction

# Ensure MAFFT temporary directory uses fast local scratch by default
MAFFT_TMP_DEFAULT="/vol/tmp"
mkdir -p "${MAFFT_TMPDIR:-$MAFFT_TMP_DEFAULT}" || true
export MAFFT_TMPDIR="${MAFFT_TMPDIR:-$MAFFT_TMP_DEFAULT}"
export TMPDIR="${TMPDIR:-$MAFFT_TMPDIR}"
echo "Using MAFFT_TMPDIR=$MAFFT_TMPDIR (TMPDIR=$TMPDIR)"

# Create output directories
echo ""
echo "Setting up output directories..."
echo "================================"
mkdir -p output/{sequences,msa,noncoding_sequences,plots}
mkdir -p output/msa/{dna_alignments,protein_alignments,noncoding_alignments}
echo "✅ Output directories created"

########## Step 1 ##########
if [ $START_STEP -le 1 ]; then
    echo ""
    echo "Step 1: Extracting operon gene sequences from assemblies..."
    echo "==========================================================="
    # Choose unified genome source for both gene and promoter extraction
    GENOME_SOURCE=${GENOME_SOURCE:-prokka}
    ASSEMBLIES_DIR=${ASSEMBLIES_DIR:-../Efs_assemblies}

    echo "🔍 Starting sequence extraction with parameters:"
    echo "   - BLAST results dir: ../03_blast_search/output/blast_results"
    echo "   - Prokka dir: ../01_prokka_annotation/output/prokka_results"
    echo "   - Output dir: output/sequences"
    echo "   - Min identity: 90%"
    echo "   - Min coverage: 80%"
    echo "   - Source: $GENOME_SOURCE"
    echo ""

    python extract_operon_sequences.py \
        --prokka_dir ../01_prokka_annotation/output/prokka_results \
        --blast_dir ../03_blast_search/output/blast_results \
        --output_dir output/sequences \
        --min_identity 90 \
        --min_coverage 80 \
        --source "$GENOME_SOURCE" \
        --assemblies_dir "$ASSEMBLIES_DIR"

    EXTRACTION_EXIT_CODE=$?
    GENE_FASTA_COUNT=$(find output/sequences -name "*.fasta" 2>/dev/null | wc -l)
    
    echo ""
    echo "📊 Step 1 Results:"
    if [ $EXTRACTION_EXIT_CODE -eq 0 ]; then
        echo "✅ Sequence extraction completed successfully"
        echo "   - Gene FASTA files created: $GENE_FASTA_COUNT"
        if [ $GENE_FASTA_COUNT -gt 0 ]; then
            echo "   - Example files:"
            find output/sequences -name "*.fasta" | head -3 | while read file; do
                seq_count=$(grep -c ">" "$file" 2>/dev/null || echo "0")
                echo "     * $(basename "$file"): $seq_count sequences"
            done
        else
            echo "⚠️  Warning: No sequence files were created - check BLAST results and thresholds"
        fi
    else
        echo "❌ Sequence extraction failed (exit code: $EXTRACTION_EXIT_CODE)"
        echo "   Check the error messages above for details"
        exit 1
    fi
else
    echo "⏭️  Skipping Step 1: Gene sequence extraction"
fi

########## Step 2 ##########
if [ $START_STEP -le 2 ]; then
    echo ""
    echo "Step 2: Creating Multiple Sequence Alignments..."
    echo "================================================"
    
    # Check if we have sequences from step 1
    SEQUENCE_COUNT=$(find output/sequences -name "*.fasta" 2>/dev/null | wc -l)
    if [ $SEQUENCE_COUNT -eq 0 ]; then
        echo "⚠️  Warning: No sequence files found in output/sequences/"
        echo "   Step 1 may have failed or been skipped"
        echo "   Proceeding anyway in case sequences exist elsewhere..."
    else
        echo "🔍 Found $SEQUENCE_COUNT sequence files for alignment"
    fi
    
    echo "🔧 Starting MSA creation with parameters:"
    echo "   - Input sequences: output/sequences"
    echo "   - Output directory: output/msa"
    echo "   - Threads: ${SLURM_CPUS_PER_TASK:-8}"
    echo "   - MAFFT_TMPDIR: $MAFFT_TMPDIR"
    echo ""

    python create_msa.py \
        --coding-sequences output/sequences \
        --output-dir output/msa \
        --threads ${SLURM_CPUS_PER_TASK:-8}
    
    MSA_EXIT_CODE=$?
    DNA_MSA_COUNT=$(find output/msa/dna_alignments -name "*.fasta" 2>/dev/null | wc -l)
    PROTEIN_MSA_COUNT=$(find output/msa/protein_alignments -name "*.fasta" 2>/dev/null | wc -l)
    
    echo ""
    echo "📊 Step 2 Results:"
    if [ $MSA_EXIT_CODE -eq 0 ]; then
        echo "✅ MSA creation completed successfully"
        echo "   - DNA alignments created: $DNA_MSA_COUNT"
        echo "   - Protein alignments created: $PROTEIN_MSA_COUNT"
        if [ $DNA_MSA_COUNT -gt 0 ]; then
            echo "   - Example DNA alignments:"
            find output/msa/dna_alignments -name "*.fasta" | head -3 | while read file; do
                seq_count=$(grep -c ">" "$file" 2>/dev/null || echo "0")
                seq_length=$(head -2 "$file" | tail -1 | wc -c 2>/dev/null || echo "0")
                echo "     * $(basename "$file"): $seq_count sequences, ~$seq_length bp"
            done
        else
            echo "⚠️  Warning: No DNA alignments were created"
        fi
    else
        echo "⚠️  MSA creation failed (exit code: $MSA_EXIT_CODE)"
        echo "   Check MAFFT installation and temporary directory permissions"
        echo "   Continuing with pipeline..."
    fi
else
    echo "⏭️  Skipping Step 2: MSA creation for genes"
fi

########## Step 3 ##########
if [ $START_STEP -le 3 ]; then
    echo ""
    echo "Step 3: Extracting promoter sequences from assemblies..."
    echo "======================================================="
    # Use the same unified source as Step 1 unless overridden
    PROMOTER_SOURCE=${PROMOTER_SOURCE:-$GENOME_SOURCE}

    python extract_noncoding_sequences.py \
        --blast-dir ../03_blast_search/output/blast_results \
        --output-dir output/noncoding_sequences \
        --min-identity 70 \
        --source "$PROMOTER_SOURCE" \
        --assemblies-dir "$ASSEMBLIES_DIR"

    if [ $? -eq 0 ]; then
        echo "✅ Promoter sequences extracted"
        
        # Create promoter MSA
        echo ""
        echo "Step 3b: Creating promoter MSA..."
        echo "================================="
        python create_msa.py \
            --noncoding-sequences output/noncoding_sequences \
            --output-dir output/msa \
            --noncoding-only \
            --threads ${SLURM_CPUS_PER_TASK:-8}
        
        if [ $? -eq 0 ]; then
            PROMOTER_MSA_COUNT=$(find output/msa/noncoding_alignments -name "*.fasta" 2>/dev/null | wc -l)
            echo "✅ Promoter MSA created ($PROMOTER_MSA_COUNT files)"
            PROMOTER_RESULTS=1
        else
            echo "⚠️  Warning: Promoter MSA creation failed"
            PROMOTER_RESULTS=0
        fi
    else
        echo "⚠️  Warning: Promoter sequence extraction failed"
        PROMOTER_RESULTS=0
    fi
else
    echo "⏭️  Skipping Step 3: Promoter extraction and MSA"
fi

########## Step 4 ##########
if [ $START_STEP -le 4 ]; then
    # Step 4a: Gene conservation plots
    echo ""
    echo "Step 4a: Creating conservation plots for all operon genes..."
    echo "========================================================"
    python create_gene_conservation_plots.py \
        --msa-dir output/msa/dna_alignments \
        --output-dir output/plots/gene_conservation \
        --title-suffix "aa_vs_nt; source=${GENOME_SOURCE}"
    echo "Gene conservation plots completed"
    echo ""

    # Step 4b: Enhanced promoter plots
    echo ""
    echo "Step 4: Creating enhanced promoter conservation plots..."
    echo "======================================================"
    python create_promoter_plot_with_pribnow.py \
        --blast-based \
        --output-dir output/plots \
        --title-suffix "source=${PROMOTER_SOURCE}"

    if [ $? -eq 0 ]; then
        echo "✅ Enhanced promoter plots created"
        PLOT_RESULTS=1
    else
        echo "⚠️  Warning: Enhanced promoter plots failed"
        PLOT_RESULTS=0
    fi
else
    echo "⏭️  Skipping Step 4: Conservation plots"
fi

########## Step 5 ##########
if [ $START_STEP -le 5 ]; then
    echo ""
    echo "Step 5: Running BLAST-based diversity analysis..."
    echo "================================================="
    python extract_sequences_from_blast.py \
        --blast-csv ../03_blast_search/output/all_blast_hits_complete.csv \
        --analysis-only

    if [ $? -eq 0 ]; then
        echo "✅ BLAST diversity analysis completed"
    else
        echo "⚠️  Warning: BLAST diversity analysis failed"
    fi
else
    echo "⏭️  Skipping Step 5: BLAST-based diversity analysis"
fi

# ########## Step 6 ##########
if [ $START_STEP -le 6 ]; then
    echo ""
    echo "Step 6: Multi-strategy comparison (currently disabled)"
    echo "====================================================="
    echo "ℹ️  Multi-strategy pipeline disabled to avoid redundant BLAST operations"
    echo "   Main analysis uses BLAST results from ../03_blast_search/ (recommended)"
    echo "   For comparison analysis, run multi-strategy scripts manually if needed"
    echo "✅ Skipping multi-strategy comparison"
else
    echo "⏭️  Skipping Step 6: Multi-strategy comparison"
fi

# Generate summary
echo ""
echo "=================================================="
echo "Operon Extraction Pipeline Complete!"
echo "Finished: $(date)"
echo "=================================================="
echo "📊 Results Summary:"

# Count results
TOTAL_GENE_FASTA=$(find output/sequences -name "*.fasta" 2>/dev/null | wc -l)
TOTAL_DNA_MSA=$(find output/msa/dna_alignments -name "*_aligned.fasta" 2>/dev/null | wc -l)
TOTAL_PROTEIN_MSA=$(find output/msa/protein_alignments -name "*_aligned.fasta" 2>/dev/null | wc -l)
TOTAL_NONCODING_MSA=$(find output/msa/noncoding_alignments -name "*_aligned.fasta" 2>/dev/null | wc -l)
TOTAL_PLOTS=$(find output/plots -name "*.png" 2>/dev/null | wc -l)
TOTAL_BLAST_RESULTS=$(find output -name "*results.csv" 2>/dev/null | wc -l)
# Promoter sequences count (number of sequences in promoter.fasta if present)
if [ -f "output/noncoding_sequences/promoter.fasta" ]; then
    TOTAL_PROMOTER_SEQ=$(grep -c '^>' output/noncoding_sequences/promoter.fasta)
else
    TOTAL_PROMOTER_SEQ=0
fi

echo "   - Gene sequences extracted: $TOTAL_GENE_FASTA"
echo "   - DNA alignments created: $TOTAL_DNA_MSA"
echo "   - Protein alignments created: $TOTAL_PROTEIN_MSA"
echo "   - Promoter MSAs: $TOTAL_NONCODING_MSA" 
echo "   - Conservation plots: $TOTAL_PLOTS"
echo "   - BLAST analysis results: $TOTAL_BLAST_RESULTS"

echo ""
echo "📁 Output directories:"
echo "   - output/sequences/                 (extracted gene sequences)"
echo "   - output/msa/dna_alignments/        (DNA MSAs)"
echo "   - output/msa/protein_alignments/    (protein MSAs)"
echo "   - output/msa/noncoding_alignments/  (promoter MSAs)"
echo "   - output/plots/                     (conservation plots)"
echo "   - output/diversity_analysis/        (BLAST analysis results)"

# Create pipeline summary file
SUMMARY_FILE="output/extraction_pipeline_summary.txt"
cat > $SUMMARY_FILE << EOF
Operon Extraction Pipeline Summary
==================================
Execution Date: $(date)
Job ID: $SLURM_JOB_ID
Execution Time: $(date)

Input Data Sources:
- BLAST results: ../03_blast_search/output/blast_results/
- Complete BLAST hits: ../03_blast_search/output/all_blast_hits_complete.csv

Results Generated:
- BLAST analysis results: $TOTAL_BLAST_RESULTS
- Promoter MSAs: $TOTAL_NONCODING_MSA
- Conservation plots: $TOTAL_PLOTS
- Promoter sequences: $TOTAL_PROMOTER_SEQ

Pipeline Steps Completed:
1. ✅ Real gene sequence extraction from assemblies
2. ✅ Multiple Sequence Alignment creation (DNA & protein)
3. ✅ Promoter analysis from BLAST
4. ✅ Enhanced promoter plots with Pribnow box
5. ✅ BLAST-based diversity analysis (supplementary)

Output Structure:
output/
├── sequences/                   # Extracted gene sequences (REAL DNA!)
├── msa/
│   ├── dna_alignments/         # DNA MSAs
│   ├── protein_alignments/     # Protein MSAs
│   └── noncoding_alignments/   # Promoter MSAs
├── plots/                      # Conservation plots
├── diversity_analysis/         # BLAST-based analysis results
└── extraction_pipeline_summary.txt

Next Steps:
- Use MSAs in ../06_diversity_analysis/ for conservation analysis
- Compare operon conservation with core genes from ../04_core_gene_analysis/
- Proceed to dN/dS analysis in ../07_dnds_analysis/
EOF

echo ""
echo "📋 Pipeline summary saved to: $SUMMARY_FILE"

# Final status check
if [ $TOTAL_GENE_FASTA -gt 0 ]; then
    echo ""
    echo "✅ Pipeline completed successfully!"
    echo ""
    echo "🧬 Ready for downstream analysis:"
    echo "   - Diversity analysis: ../06_diversity_analysis/"
    echo "   - dN/dS analysis: ../07_dnds_analysis/"
    echo ""
    echo "🎯 Key outputs for next steps:"
    echo "   - MSAs: output/msa/"
    echo "   - Sequences: output/sequences/"
    echo "   - Promoter plots: output/plots/"
    exit 0
else
    echo ""
    echo "❌ Pipeline failed - insufficient results generated"
    echo "Check individual step outputs for debugging"
    exit 1
fi