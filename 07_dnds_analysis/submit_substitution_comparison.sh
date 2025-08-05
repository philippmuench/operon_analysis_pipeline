#!/bin/bash
#SBATCH --job-name=substitution_comp
#SBATCH --output=logs/substitution_comparison_%j.out
#SBATCH --error=logs/substitution_comparison_%j.err
#SBATCH --time=06:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --partition=cpu

echo "=== Substitution Analysis: Operon vs Core Genes ==="
echo "Started at: $(date)"
echo

# Load conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate efs_diversity

# Run the comparison
echo "Analyzing substitution patterns in operon and core genes..."
python compare_substitutions_all.py

# Check if successful
if [ $? -eq 0 ]; then
    echo -e "\nAnalysis completed successfully at $(date)"
    
    # Show summary
    if [ -f "substitution_comparison_results.csv" ]; then
        echo -e "\nQuick summary:"
        python -c "
import pandas as pd
df = pd.read_csv('substitution_comparison_results.csv')
operon = df[df['gene_type'] == 'operon']
core = df[df['gene_type'] == 'core']

print(f'Genes analyzed:')
print(f'  Operon: {len(operon)}')
print(f'  Core: {len(core)}')
print(f'\\nMean dN/dS ratios:')
print(f'  Operon: {operon[\"dn_ds_ratio\"][operon[\"dn_ds_ratio\"] < 10].mean():.3f}')
print(f'  Core: {core[\"dn_ds_ratio\"][core[\"dn_ds_ratio\"] < 10].mean():.3f}')
"
    fi
    
    echo -e "\nOutput files:"
    ls -lh operon_vs_core_substitutions.png substitution_comparison_results.csv 2>/dev/null
else
    echo "Error: Analysis failed"
    exit 1
fi