#!/usr/bin/env python3
"""
Create a comprehensive summary of the core gene analysis results.
"""

import pandas as pd
import os

def main():
    print("📊 CORE GENE ANALYSIS SUMMARY")
    print("==============================")
    print()
    
    # Check for output files
    output_dir = "output"
    
    # Gene prevalence stats
    prevalence_file = os.path.join(output_dir, "gene_prevalence_stats.csv")
    if os.path.exists(prevalence_file):
        prevalence_df = pd.read_csv(prevalence_file)
        print(f"✅ GENE PREVALENCE ANALYSIS:")
        print(f"   📊 Total unique genes analyzed: {len(prevalence_df):,}")
        print(f"   📈 Prevalence range: {prevalence_df['prevalence'].min():.2f} - {prevalence_df['prevalence'].max():.2f}")
        print()
    
    # Threshold analysis
    threshold_file = os.path.join(output_dir, "core_gene_threshold_summary.csv")
    if os.path.exists(threshold_file):
        threshold_df = pd.read_csv(threshold_file)
        print(f"✅ THRESHOLD ANALYSIS:")
        print(f"   📊 Core genes at 100% threshold: {threshold_df[threshold_df['Threshold (%)'] == 100]['Core genes'].iloc[0]:,}")
        print(f"   📊 Core genes at 95% threshold: {threshold_df[threshold_df['Threshold (%)'] == 95]['Core genes'].iloc[0]:,}")
        print(f"   📊 Core genes at 90% threshold: {threshold_df[threshold_df['Threshold (%)'] == 90]['Core genes'].iloc[0]:,}")
        print()
    
    # Conservation metrics
    conservation_file = os.path.join(output_dir, "core_gene_conservation_summary.txt")
    if os.path.exists(conservation_file):
        print(f"✅ CONSERVATION ANALYSIS:")
        with open(conservation_file, 'r') as f:
            lines = f.readlines()
            for line in lines[:10]:  # Show first 10 lines
                print(f"   {line.strip()}")
        print()
    
    # Generated plots
    plot_files = [
        ("Core gene threshold curve", "core_gene_threshold_curve.png"),
        ("Core gene threshold curve (zoomed)", "core_gene_threshold_curve_zoomed.png")
    ]
    
    plots_found = []
    for name, filename in plot_files:
        filepath = os.path.join(output_dir, filename)
        if os.path.exists(filepath):
            plots_found.append(f"   📊 {filepath}")
    
    if plots_found:
        print(f"✅ GENERATED PLOTS:")
        for plot in plots_found:
            print(plot)
        print()
    
    print("🎯 ANALYSIS COMPLETE!")
    print("Ready for scientific interpretation and publication.")

if __name__ == "__main__":
    main()
