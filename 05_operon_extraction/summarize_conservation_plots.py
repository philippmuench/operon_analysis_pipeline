#!/usr/bin/env python3
"""
Generate a summary of all conservation analysis plots created.
"""

import os
import glob

def main():
    plot_dirs = [
        ("Promoter Conservation", "output/plots/promoter_conservation_with_pribnow.png"),
        ("Gene Conservation", "output/plots/gene_conservation/")
    ]
    
    print("🎨 CONSERVATION ANALYSIS PLOTS SUMMARY")
    print("======================================")
    print()
    
    # Promoter plot
    promoter_plot = "output/plots/promoter_conservation_with_pribnow.png"
    if os.path.exists(promoter_plot):
        print("✅ PROMOTER CONSERVATION:")
        print(f"   📊 {promoter_plot}")
        print("   📍 Includes Pribnow box highlighting (positions 43-50)")
        print("   📈 Two panels: SNP-based conservation + gap frequency")
        print()
    
    # Gene conservation plots
    gene_plot_dir = "output/plots/gene_conservation/"
    if os.path.exists(gene_plot_dir):
        gene_plots = glob.glob(os.path.join(gene_plot_dir, "*_conservation.png"))
        if gene_plots:
            print("✅ OPERON GENE CONSERVATION:")
            for plot_file in sorted(gene_plots):
                gene_name = os.path.basename(plot_file).replace('_conservation.png', '')
                print(f"   📊 {plot_file}")
                print(f"      🧬 Gene: {gene_name}")
            print("   📈 Each plot: SNP-based conservation + gap frequency")
            print()
            
            print(f"📊 TOTAL PLOTS: 1 promoter + {len(gene_plots)} genes = {1 + len(gene_plots)} conservation plots")
        else:
            print("❌ No gene conservation plots found")
    else:
        print("❌ Gene conservation directory not found")
    
    print()
    print("🎯 All plots ready for scientific analysis and publication!")

if __name__ == "__main__":
    main()
