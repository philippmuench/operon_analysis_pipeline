#!/bin/bash

# Check Prokka annotation progress

OUTPUT_DIR="../prokka_output"
GENOME_DIR="../Efs_assemblies"

# Count total genomes
TOTAL_GENOMES=$(ls $GENOME_DIR/*.fasta.gz 2>/dev/null | wc -l)

# Count completed genomes (those with .gff files)
COMPLETED=$(find $OUTPUT_DIR -name "*.gff" 2>/dev/null | wc -l)

echo "Prokka Annotation Progress"
echo "========================="
echo "Total genomes: $TOTAL_GENOMES"
echo "Completed: $COMPLETED"
echo "Remaining: $(($TOTAL_GENOMES - $COMPLETED))"
echo "Progress: $(awk "BEGIN {printf \"%.1f\", $COMPLETED/$TOTAL_GENOMES*100}")%"