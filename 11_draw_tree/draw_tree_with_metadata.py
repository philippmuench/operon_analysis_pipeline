#!/usr/bin/env python3
"""
Complete pipeline to create GraPhLan tree visualization with metadata coloring.
Reads tree file and metadata, creates annotation, and generates visualization.
"""

import os
import sys
import re
import pandas as pd
import subprocess
from collections import Counter

def extract_samples_from_tree(tree_file):
    """Extract sample names from Newick tree file."""
    print("Reading tree file...")
    with open(tree_file, 'r') as f:
        tree_content = f.read()
    
    # Extract sample names (format: ENT_XXXXXXAA_AS.result)
    pattern = r'ENT_[A-Z]+\d+[A-Z]+_AS\.result'
    samples = list(set(re.findall(pattern, tree_content)))
    print(f"  Found {len(samples)} unique samples in tree")
    return samples

def read_metadata(metadata_file):
    """Read metadata file and create mappings."""
    print("Reading metadata...")
    df = pd.read_csv(metadata_file, sep='\t')
    print(f"  Loaded {len(df)} rows from metadata")
    
    # Create AS_barcode to niche mapping
    barcode_to_niche = {}
    for _, row in df.iterrows():
        if pd.notna(row['AS_barcode']) and pd.notna(row['Source Niche']):
            barcode_to_niche[row['AS_barcode']] = row['Source Niche']
        elif pd.notna(row['AS_barcode']):
            barcode_to_niche[row['AS_barcode']] = 'Unknown'
    
    return barcode_to_niche

def create_color_palette():
    """Create scientific color palette for niches."""
    return {
        'Human': '#E74C3C',           # Red
        'Livestock': '#27AE60',        # Green  
        'Poultry': '#3498DB',          # Blue
        'Wild Animal': '#9B59B6',      # Purple
        'Wild animal': '#9B59B6',      # Purple (alternative capitalization)
        'Companion Animal': '#F39C12', # Orange
        'Environment': '#1ABC9C',      # Turquoise
        'Food': '#E67E22',             # Dark Orange
        'Animal Feed': '#16A085',      # Dark Turquoise
        'Laboratory': '#8E44AD',       # Dark Purple
        'Aquatic Animal': '#2980B9',   # Dark Blue
        'Unknown': '#95A5A6',          # Gray
        'ND': '#BDC3C7'                # Light Gray
    }

def create_graphlan_annotation(tree_samples, metadata, output_file):
    """Create GraPhLan annotation file."""
    print("Creating annotation file...")
    
    colors = create_color_palette()
    matched = 0
    niche_counts = Counter()
    
    lines = []
    
    # Header and global settings
    lines.extend([
        "# GraPhLan annotation for E. faecalis phylogeny",
        "title\tE. faecalis phylogeny by niche",
        "title_font_size\t25",
        "total_plotted_degrees\t340",
        "start_rotation\t90",
        "branch_thickness\t1.0",
        "branch_color\t#606060",
        "clade_marker_size\t3.0",
        "clade_marker_edge_width\t0.5",
        "clade_marker_edge_color\t#303030",
        ""
    ])
    
    # Add sample annotations
    for sample in tree_samples:
        # Extract barcode from sample name
        barcode = sample.replace('.result', '') if sample.endswith('.result') else sample
        
        if barcode in metadata:
            niche = metadata[barcode]
            color = colors.get(niche, '#95A5A6')
            lines.append(f"{sample}\tclade_marker_color\t{color}")
            niche_counts[niche] += 1
            matched += 1
    
    # Write file
    with open(output_file, 'w') as f:
        f.write('\n'.join(lines))
    
    print(f"  Matched {matched}/{len(tree_samples)} samples")
    print(f"  Created annotation with {len(lines)} lines")
    
    # Print summary
    print("\nNiche distribution:")
    total = sum(niche_counts.values())
    for niche, count in niche_counts.most_common(10):
        pct = (count/total)*100
        print(f"  {niche:20s}: {count:5d} ({pct:5.1f}%)")
    
    return output_file

def run_graphlan(tree_file, annotation_file, output_prefix):
    """Run GraPhLan to create visualization."""
    print("\nRunning GraPhLan visualization...")
    
    # Step 1: Create annotated XML
    xml_file = f"{output_prefix}.xml"
    cmd1 = ["graphlan_annotate.py", tree_file, xml_file, "--annot", annotation_file]
    
    print("  Creating annotated XML...")
    result = subprocess.run(cmd1, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: Failed to create XML\n{result.stderr}")
        return False
    
    # Step 2: Generate PNG visualization
    png_file = f"{output_prefix}.png"
    cmd2 = ["graphlan.py", xml_file, png_file, 
            "--size", "20", "--dpi", "300", "--pad", "0.2"]
    
    print("  Generating PNG visualization...")
    result = subprocess.run(cmd2, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: Failed to create PNG\n{result.stderr}")
        return False
    
    print(f"  Successfully created: {png_file}")
    
    # Step 3: Generate SVG for publication (optional)
    svg_file = f"{output_prefix}.svg"
    cmd3 = ["graphlan.py", xml_file, svg_file,
            "--size", "20", "--dpi", "300", "--pad", "0.2", "--format", "svg"]
    
    print("  Generating SVG for publication...")
    result = subprocess.run(cmd3, capture_output=True, text=True)
    if result.returncode == 0:
        print(f"  Successfully created: {svg_file}")
    else:
        print("  Warning: SVG creation failed, but PNG is available")
    
    return True

def main():
    """Main pipeline function."""
    
    # Configuration
    base_dir = "/vol/projects/BIFO/genomenet/baerbel_science_rebuttal/operon_analysis"
    tree_file = os.path.join(base_dir, "10_phylophlan/output_isolates/RAxML_bestTree.input_isolates_refined.tre")
    metadata_file = os.path.join(base_dir, "00_annotation/8587_Efs_metadata_ASbarcode.txt")
    work_dir = os.path.join(base_dir, "11_draw_tree")
    output_dir = os.path.join(work_dir, "output")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Change to working directory
    os.chdir(work_dir)
    
    print("="*60)
    print("GraPhLan Tree Visualization Pipeline")
    print("="*60)
    print()
    
    # Check input files
    if not os.path.exists(tree_file):
        print(f"ERROR: Tree file not found: {tree_file}")
        return 1
    if not os.path.exists(metadata_file):
        print(f"ERROR: Metadata file not found: {metadata_file}")
        return 1
    
    # Run pipeline
    try:
        # Extract samples from tree
        tree_samples = extract_samples_from_tree(tree_file)
        
        # Read metadata
        metadata = read_metadata(metadata_file)
        
        # Create annotation file in output directory
        annotation_file = os.path.join(output_dir, "graphlan_annotation.txt")
        create_graphlan_annotation(tree_samples, metadata, annotation_file)
        
        # Run GraPhLan with output in output directory
        output_prefix = os.path.join(output_dir, "Efs_phylogeny_niche")
        success = run_graphlan(tree_file, annotation_file, output_prefix)
        
        if success:
            print("\n" + "="*60)
            print("SUCCESS! Visualization complete")
            print("="*60)
            print(f"\nOutput files in output/ directory:")
            print(f"  Efs_phylogeny_niche.png - High-resolution tree image")
            print(f"  Efs_phylogeny_niche.svg - Vector format for publication")
            print(f"  Efs_phylogeny_niche.xml - Annotated tree data")
            return 0
        else:
            print("\nERROR: Visualization failed")
            return 1
            
    except Exception as e:
        print(f"\nERROR: Pipeline failed with exception: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())