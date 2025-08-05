#!/bin/bash

# Validate Prokka annotation outputs
# Checks that all expected files exist and are non-empty

# Parse command line arguments
TEST_MODE=false
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --test) TEST_MODE=true ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Set paths based on mode
if [ "$TEST_MODE" = true ]; then
    OUTPUT_DIR="output/prokka_results"
    GENOME_LIST_FILE="genome_list_test.txt"
    LOG_FILE="output/validation_report_test_$(date +%Y%m%d_%H%M%S).log"
else
    OUTPUT_DIR="../prokka_output"
    GENOME_LIST_FILE="genome_list.txt"
    LOG_FILE="validation_report_$(date +%Y%m%d_%H%M%S).log"
fi

GENOME_DIR="../Efs_assemblies"

# Expected Prokka output files per genome
EXPECTED_FILES=("gff" "faa" "ffn" "fna" "gbk" "log" "txt" "tsv")

# Initialize counters
TOTAL_GENOMES=0
COMPLETE_GENOMES=0
INCOMPLETE_GENOMES=0
MISSING_OUTPUTS=0
EMPTY_FILES=0

echo "Prokka Output Validation Report" | tee $LOG_FILE
echo "===============================" | tee -a $LOG_FILE
echo "Started: $(date)" | tee -a $LOG_FILE
echo "" | tee -a $LOG_FILE

# Get list of genomes based on mode
if [ "$TEST_MODE" = true ]; then
    # Read from genome_list_test.txt file
    if [ -f "$GENOME_LIST_FILE" ]; then
        GENOME_LIST=$(cat "$GENOME_LIST_FILE" | sed 's|^.*/||' | sed 's|^|'$GENOME_DIR'/|')
        TOTAL_GENOMES=$(cat "$GENOME_LIST_FILE" | wc -l)
    else
        echo "Error: Test genome list file not found: $GENOME_LIST_FILE"
        exit 1
    fi
else
    GENOME_LIST=$(ls $GENOME_DIR/*.fasta.gz 2>/dev/null)
    TOTAL_GENOMES=$(echo "$GENOME_LIST" | wc -l)
fi

echo "Total genomes to process: $TOTAL_GENOMES" | tee -a $LOG_FILE
echo "" | tee -a $LOG_FILE

# Check each genome
echo "Checking individual genomes..." | tee -a $LOG_FILE
echo "------------------------------" | tee -a $LOG_FILE

for GENOME in $GENOME_LIST; do
    BASENAME=$(basename "$GENOME" .fasta.gz)
    GENOME_DIR="$OUTPUT_DIR/$BASENAME"
    GENOME_COMPLETE=true
    GENOME_ISSUES=""
    
    # Check if output directory exists
    if [ ! -d "$GENOME_DIR" ]; then
        echo "ERROR: Missing output directory for $BASENAME" | tee -a $LOG_FILE
        INCOMPLETE_GENOMES=$((INCOMPLETE_GENOMES + 1))
        MISSING_OUTPUTS=$((MISSING_OUTPUTS + ${#EXPECTED_FILES[@]}))
        continue
    fi
    
    # Check each expected file
    for EXT in "${EXPECTED_FILES[@]}"; do
        FILE="$GENOME_DIR/$BASENAME.$EXT"
        
        if [ ! -f "$FILE" ]; then
            GENOME_COMPLETE=false
            GENOME_ISSUES="${GENOME_ISSUES}  - Missing: $BASENAME.$EXT\n"
            MISSING_OUTPUTS=$((MISSING_OUTPUTS + 1))
        elif [ ! -s "$FILE" ]; then
            GENOME_COMPLETE=false
            GENOME_ISSUES="${GENOME_ISSUES}  - Empty: $BASENAME.$EXT\n"
            EMPTY_FILES=$((EMPTY_FILES + 1))
        fi
    done
    
    # Report genome status
    if [ "$GENOME_COMPLETE" = true ]; then
        COMPLETE_GENOMES=$((COMPLETE_GENOMES + 1))
    else
        INCOMPLETE_GENOMES=$((INCOMPLETE_GENOMES + 1))
        echo "INCOMPLETE: $BASENAME" | tee -a $LOG_FILE
        echo -e "$GENOME_ISSUES" | tee -a $LOG_FILE
    fi
done

# Summary statistics
echo "" | tee -a $LOG_FILE
echo "Summary Statistics" | tee -a $LOG_FILE
echo "==================" | tee -a $LOG_FILE
echo "Total genomes: $TOTAL_GENOMES" | tee -a $LOG_FILE
echo "Complete genomes: $COMPLETE_GENOMES ($(awk "BEGIN {printf \"%.1f\", $COMPLETE_GENOMES/$TOTAL_GENOMES*100}")%)" | tee -a $LOG_FILE
echo "Incomplete genomes: $INCOMPLETE_GENOMES" | tee -a $LOG_FILE
echo "Missing files: $MISSING_OUTPUTS" | tee -a $LOG_FILE
echo "Empty files: $EMPTY_FILES" | tee -a $LOG_FILE
echo "" | tee -a $LOG_FILE

# Additional checks for complete genomes
if [ $COMPLETE_GENOMES -gt 0 ]; then
    echo "Additional Statistics for Complete Genomes" | tee -a $LOG_FILE
    echo "==========================================" | tee -a $LOG_FILE
    
    # Average file sizes
    echo "Average file sizes:" | tee -a $LOG_FILE
    for EXT in "${EXPECTED_FILES[@]}"; do
        AVG_SIZE=$(find $OUTPUT_DIR -name "*.$EXT" -exec stat -c%s {} \; | \
                   awk '{sum+=$1; count++} END {if(count>0) printf "%.0f", sum/count; else print "0"}')
        AVG_SIZE_MB=$(awk "BEGIN {printf \"%.2f\", $AVG_SIZE/1048576}")
        echo "  .$EXT files: ${AVG_SIZE_MB} MB" | tee -a $LOG_FILE
    done
    
    # Check for potential issues in GFF files
    echo "" | tee -a $LOG_FILE
    echo "Checking GFF file contents..." | tee -a $LOG_FILE
    
    # Count genomes with at least one CDS
    GENOMES_WITH_CDS=$(find $OUTPUT_DIR -name "*.gff" -exec grep -l "^[^#].*CDS" {} \; | wc -l)
    echo "  Genomes with CDS annotations: $GENOMES_WITH_CDS" | tee -a $LOG_FILE
    
    # Sample CDS counts from a few genomes
    echo "  Sample CDS counts:" | tee -a $LOG_FILE
    find $OUTPUT_DIR -name "*.gff" | head -5 | while read GFF; do
        BASENAME=$(basename "$GFF" .gff)
        CDS_COUNT=$(grep -c "^[^#].*CDS" "$GFF" || echo "0")
        echo "    $BASENAME: $CDS_COUNT CDS features" | tee -a $LOG_FILE
    done
fi

# Generate list of incomplete genomes for easy re-processing
if [ $INCOMPLETE_GENOMES -gt 0 ]; then
    INCOMPLETE_LIST="incomplete_genomes_$(date +%Y%m%d_%H%M%S).txt"
    echo "" | tee -a $LOG_FILE
    echo "Generating list of incomplete genomes: $INCOMPLETE_LIST" | tee -a $LOG_FILE
    
    > $INCOMPLETE_LIST
    for GENOME in $GENOME_LIST; do
        BASENAME=$(basename "$GENOME" .fasta.gz)
        GENOME_DIR="$OUTPUT_DIR/$BASENAME"
        
        # Check if this genome is incomplete
        if [ ! -d "$GENOME_DIR" ]; then
            echo "$GENOME" >> $INCOMPLETE_LIST
            continue
        fi
        
        INCOMPLETE=false
        for EXT in "${EXPECTED_FILES[@]}"; do
            FILE="$GENOME_DIR/$BASENAME.$EXT"
            if [ ! -f "$FILE" ] || [ ! -s "$FILE" ]; then
                INCOMPLETE=true
                break
            fi
        done
        
        if [ "$INCOMPLETE" = true ]; then
            echo "$GENOME" >> $INCOMPLETE_LIST
        fi
    done
    
    echo "Total incomplete genomes in list: $(wc -l < $INCOMPLETE_LIST)" | tee -a $LOG_FILE
fi

echo "" | tee -a $LOG_FILE
echo "Validation completed: $(date)" | tee -a $LOG_FILE
echo "Full report saved to: $LOG_FILE" | tee -a $LOG_FILE

# Exit with error if any genomes are incomplete
if [ $INCOMPLETE_GENOMES -gt 0 ]; then
    exit 1
else
    echo "All genomes successfully processed!" | tee -a $LOG_FILE
    exit 0
fi