# Tree Visualization with GraPhLan

Visualizes the E. faecalis phylogenetic tree colored by source niche metadata.

## Usage

```bash
# Submit SLURM job
sbatch submit_tree_visualization.sh

# Or run directly
python3 draw_tree_with_metadata.py
```

## Input
- Tree: `../10_phylophlan/output_isolates/RAxML_bestTree.input_isolates_refined.tre`
- Metadata: `../00_annotation/8587_Efs_metadata_ASbarcode.txt`

## Output
Files generated in `output/` directory:
- `Efs_phylogeny_niche.png` - Tree visualization (300 DPI)
- `Efs_phylogeny_niche.svg` - Vector format for publication
- `graphlan_annotation.txt` - GraPhLan annotation file

## Requirements
- GraPhLan (in `efs_diversity` conda environment)
- Python with pandas