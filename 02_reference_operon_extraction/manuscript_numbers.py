#!/usr/bin/env python3
"""
Generate manuscript statistics for reference operon extraction methods section.
"""

import os
from Bio import SeqIO

def analyze_genbank_file():
    """Analyze the input GenBank file."""
    gb_file = "operon.gb"
    if not os.path.exists(gb_file):
        return {"error": "GenBank file not found"}
    
    stats = {}
    
    try:
        with open(gb_file, 'r') as f:
            for record in SeqIO.parse(f, 'genbank'):
                stats['accession'] = record.id
                stats['description'] = record.description
                stats['sequence_length'] = len(record.seq)
                stats['organism'] = record.annotations.get('organism', 'Unknown')
                
                # Count features
                feature_counts = {}
                for feature in record.features:
                    ftype = feature.type
                    feature_counts[ftype] = feature_counts.get(ftype, 0) + 1
                
                stats['features'] = feature_counts
                break  # Only process first record
    except Exception as e:
        stats['error'] = f"Error reading GenBank file: {e}"
    
    return stats

def analyze_output_files():
    """Analyze the generated output files."""
    output_dir = "output"
    stats = {}
    
    # Check output files
    output_files = [
        ("operon_genes.tsv", "Gene information table"),
        ("operon_genes_protein.fasta", "Protein sequences"),
        ("operon_genes_nt.fasta", "Nucleotide gene sequences"),
        ("operon_noncoding_nt.fasta", "Non-coding sequences")
    ]
    
    for filename, description in output_files:
        filepath = os.path.join(output_dir, filename)
        if os.path.exists(filepath):
            stats[filename] = {"exists": True, "description": description}
            
            # Count sequences in FASTA files
            if filename.endswith('.fasta'):
                try:
                    with open(filepath, 'r') as f:
                        seq_count = sum(1 for line in f if line.startswith('>'))
                    stats[filename]['sequence_count'] = seq_count
                    
                    # Get sequence lengths
                    lengths = []
                    with open(filepath, 'r') as f:
                        for record in SeqIO.parse(f, 'fasta'):
                            lengths.append(len(record.seq))
                    
                    if lengths:
                        stats[filename]['avg_length'] = sum(lengths) / len(lengths)
                        stats[filename]['min_length'] = min(lengths)
                        stats[filename]['max_length'] = max(lengths)
                        stats[filename]['total_length'] = sum(lengths)
                
                except Exception as e:
                    stats[filename]['error'] = str(e)
            
            # Count lines in TSV
            elif filename.endswith('.tsv'):
                try:
                    with open(filepath, 'r') as f:
                        line_count = sum(1 for line in f) - 1  # Subtract header
                    stats[filename]['gene_count'] = line_count
                except Exception as e:
                    stats[filename]['error'] = str(e)
        else:
            stats[filename] = {"exists": False}
    
    return stats

def analyze_gene_content():
    """Analyze the specific genes in the operon."""
    gene_info = {
        'operon_name': 'Fructoselysine/glucoselysine utilization operon',
        'genes': [
            {'name': 'frpC', 'product': 'fructoselysine-6-phosphate deglycase', 'type': 'CDS'},
            {'name': 'glpC', 'product': 'glucoselysine-6-phosphate deglycase', 'type': 'CDS'},
            {'name': 'ptsD', 'product': 'PTS system mannose/fructose/sorbose family IID component', 'type': 'CDS'},
            {'name': 'ptsC', 'product': 'PTS system sorbose-specific IIC component', 'type': 'CDS'},
            {'name': 'ptsB', 'product': 'PTS system sorbose subfamily IIB component', 'type': 'CDS'},
            {'name': 'ptsA', 'product': 'PTS system fructose IIA component', 'type': 'CDS'},
            {'name': 'fruR', 'product': 'sigma-54 dependent transcriptional regulator', 'type': 'CDS'},
            {'name': 'promoter', 'product': 'operon promoter region', 'type': 'regulatory'},
            {'name': 'pribnow_box', 'product': 'Pribnow box (-10 box)', 'type': 'regulatory'}
        ],
        'total_genes': 7,
        'total_regulatory': 2,
        'pathway': 'Fructoselysine and glucoselysine utilization via PTS transport system'
    }
    
    return gene_info

def generate_manuscript_stats(output_file=None):
    """Generate all statistics for manuscript.
    
    Args:
        output_file (str): Optional path to save output to file. If None, prints to console.
    """
    
    # Collect all output lines
    output_lines = []
    
    def add_line(text=""):
        """Add a line to both console and file output."""
        print(text)
        output_lines.append(text)
    
    add_line("Reference Operon Extraction Statistics for Manuscript")
    add_line("=" * 60)
    
    # GenBank input analysis
    add_line("\n1. Reference Genome Source:")
    gb_stats = analyze_genbank_file()
    if "error" not in gb_stats:
        add_line(f"   Organism: {gb_stats.get('organism', 'Unknown')}")
        add_line(f"   Accession: {gb_stats.get('accession', 'Unknown')}")
        add_line(f"   Sequence length: {gb_stats.get('sequence_length', 0):,} bp")
        if 'features' in gb_stats:
            add_line(f"   GenBank features: {sum(gb_stats['features'].values())} total")
            for ftype, count in sorted(gb_stats['features'].items()):
                add_line(f"     {ftype}: {count}")
    else:
        add_line(f"   {gb_stats['error']}")
    
    # Gene content analysis
    add_line("\n2. Operon Gene Content:")
    gene_info = analyze_gene_content()
    add_line(f"   Operon: {gene_info['operon_name']}")
    add_line(f"   Coding genes: {gene_info['total_genes']}")
    add_line(f"   Regulatory elements: {gene_info['total_regulatory']}")
    add_line(f"   Biological pathway: {gene_info['pathway']}")
    
    add_line("\n   Gene composition:")
    for gene in gene_info['genes']:
        if gene['type'] == 'CDS':
            add_line(f"     {gene['name']}: {gene['product']}")
    
    add_line("\n   Regulatory elements:")
    for gene in gene_info['genes']:
        if gene['type'] == 'regulatory':
            add_line(f"     {gene['name']}: {gene['product']}")
    
    # Output files analysis
    add_line("\n3. Generated Query Sequences:")
    output_stats = analyze_output_files()
    
    for filename, stats in output_stats.items():
        if stats['exists']:
            add_line(f"   {filename}:")
            add_line(f"     Description: {stats['description']}")
            
            if 'sequence_count' in stats:
                add_line(f"     Sequences: {stats['sequence_count']}")
                add_line(f"     Average length: {stats['avg_length']:.0f} residues/bp")
                add_line(f"     Length range: {stats['min_length']}-{stats['max_length']}")
                add_line(f"     Total length: {stats['total_length']:,} residues/bp")
            
            if 'gene_count' in stats:
                add_line(f"     Genes catalogued: {stats['gene_count']}")
        else:
            add_line(f"   {filename}: Not found")
    
    # Summary for manuscript
    add_line("\n" + "=" * 60)
    add_line("MANUSCRIPT NUMBERS SUMMARY:")
    add_line("=" * 60)
    
    protein_stats = output_stats.get('operon_genes_protein.fasta', {})
    nt_stats = output_stats.get('operon_genes_nt.fasta', {})
    noncoding_stats = output_stats.get('operon_noncoding_nt.fasta', {})
    
    add_line(f"Reference organism: {gb_stats.get('organism', 'E. faecalis V583')}")
    add_line(f"Operon genes extracted: {gene_info['total_genes']}")
    add_line(f"Regulatory elements: {gene_info['total_regulatory']}")
    
    if protein_stats.get('sequence_count'):
        add_line(f"Protein query sequences: {protein_stats['sequence_count']}")
        add_line(f"Average protein length: {protein_stats['avg_length']:.0f} amino acids")
    
    if nt_stats.get('sequence_count'):
        add_line(f"Nucleotide gene sequences: {nt_stats['sequence_count']}")
        add_line(f"Average gene length: {nt_stats['avg_length']:.0f} bp")
    
    if noncoding_stats.get('sequence_count'):
        add_line(f"Regulatory sequences: {noncoding_stats['sequence_count']}")
    
    add_line(f"Output formats: TSV table + 3 FASTA files (protein, nucleotide, regulatory)")
    
    # Write to file if specified
    if output_file:
        try:
            with open(output_file, 'w') as f:
                for line in output_lines:
                    f.write(line + '\n')
            add_line(f"\nResults saved to: {output_file}")
        except Exception as e:
            add_line(f"\nError saving to file: {e}")
    
    return {
        'genbank_stats': gb_stats,
        'gene_info': gene_info,
        'output_stats': output_stats,
        'output_lines': output_lines
    }

if __name__ == "__main__":
    import sys
    
    # Change to script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    
    # Check for output file argument
    output_file = None
    if len(sys.argv) > 1:
        output_file = sys.argv[1]
    
    stats = generate_manuscript_stats(output_file)
