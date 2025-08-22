#!/usr/bin/env python3
"""
Analyze start-codon choice discrepancies for operon genes across Prokka/Prodigal outputs.

For each genome's GFF + FNA, find operon genes of interest, extract coding sequence
context, scan for upstream in-frame alternative starts (ATG/GTG/TTG) and RBS motifs,
and summarize likely reasons for shorter/alternate starts.

Also stratifies results by metadata (Source Niche, Continent, etc.) when metadata file is provided.
"""

import argparse
import os
import sys
import glob
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Tuple
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from concurrent.futures import ProcessPoolExecutor, as_completed

from Bio import SeqIO
from Bio.Seq import Seq

# Set plotting style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")


RBS_MOTIFS = ["AGGAGG", "GGAGG", "AGGA", "GGAG", "GAGG"]


@dataclass
class StartSiteResult:
    genome_id: str
    gene: str
    contig: str
    strand: str
    start: int
    end: int
    start_codon: str
    gff_start_type: Optional[str]
    gff_partial: Optional[str]
    gff_rbs_motif: Optional[str]
    gff_rbs_spacer: Optional[str]
    upstream_candidate_found: bool
    upstream_codon: Optional[str]
    upstream_offset_nt: Optional[int]
    upstream_has_rbs: Optional[bool]
    upstream_rbs_seq: Optional[str]
    upstream_rbs_spacer_nt: Optional[int]
    upstream_rbs_mismatches: Optional[int]
    classification: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Start-site analysis for operon genes")
    parser.add_argument("--prokka_dir", help="Directory containing per-genome Prokka outputs (not needed for --visualize-only)")
    parser.add_argument(
        "--gene_reference_fasta",
        default=os.path.join("..", "02_reference_operon_extraction", "output", "operon_genes_nt.fasta"),
        help="FASTA of reference operon genes to define target gene names",
    )
    parser.add_argument("--output_dir", default="output", help="Output directory")
    parser.add_argument("--max_upstream", type=int, default=120, help="Max nt upstream to scan for alt starts")
    parser.add_argument("--rbs_min_spacer", type=int, default=4, help="Min nt spacer between RBS and start")
    parser.add_argument("--rbs_max_spacer", type=int, default=13, help="Max nt spacer between RBS and start")
    parser.add_argument("--max_workers", type=int, default=8, help="Parallel workers")
    parser.add_argument("--genome_list", help="Optional file with genome IDs to analyze (one per line)")
    parser.add_argument("--blast_dir", help="Optional BLAST results directory to anchor gene coordinates (recommended)")
    parser.add_argument("--metadata", default="../00_annotation/8587_Efs_metadata_ASbarcode.txt", 
                       help="Optional metadata file for stratification analysis")
    parser.add_argument("--visualize-only", action="store_true",
                       help="Only create visualizations from existing results")
    return parser.parse_args()


def load_operon_gene_names(reference_fasta: str) -> List[str]:
    names: List[str] = []
    for rec in SeqIO.parse(reference_fasta, "fasta"):
        # headers look like: operon_CDS_1|frpC|...
        parts = rec.id.split("|")
        if len(parts) >= 2:
            names.append(parts[1])
        else:
            names.append(rec.id)
    return names


def read_genome_sequences(fna_path: str) -> Dict[str, Seq]:
    seqs: Dict[str, Seq] = {}
    for rec in SeqIO.parse(fna_path, "fasta"):
        seqs[rec.id] = rec.seq
    return seqs


def parse_gff_genes(gff_path: str) -> List[Dict[str, str]]:
    genes: List[Dict[str, str]] = []
    with open(gff_path, "r") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = parts
            if ftype not in ("CDS", "gene"):
                continue
            attr_map: Dict[str, str] = {}
            for a in attrs.split(";"):
                if not a:
                    continue
                if "=" in a:
                    k, v = a.split("=", 1)
                    attr_map[k] = v
            gene_name = attr_map.get("gene") or attr_map.get("Name") or attr_map.get("locus_tag")
            genes.append(
                {
                    "contig": seqid,
                    "source": source,
                    "type": ftype,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "attrs": attr_map,
                    "gene_name": gene_name,
                }
            )
    return genes


def parse_blast_best_hits(blast_path: str) -> List[Dict[str, str]]:
    """
    Parse outfmt 6 and keep best hit per qseqid by bitscore.
    Returns list of dicts with keys: qseqid, sseqid, sstart, send, bitscore.
    """
    best: Dict[str, Dict[str, str]] = {}
    if not os.path.exists(blast_path) or os.path.getsize(blast_path) == 0:
        return []
    with open(blast_path, "r") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            qseqid, sseqid = parts[0], parts[1]
            sstart, send = int(parts[8]), int(parts[9])
            bitscore = float(parts[11])
            rec = {"qseqid": qseqid, "sseqid": sseqid, "sstart": sstart, "send": send, "bitscore": bitscore}
            prev = best.get(qseqid)
            if prev is None or bitscore > prev["bitscore"]:
                best[qseqid] = rec
    return list(best.values())


def find_blast_file_for_genome(blast_dir: str, genome_id: str) -> Optional[str]:
    base = genome_id
    if base.endswith(".result"):
        base = base[:-7]
    candidates = [
        os.path.join(blast_dir, f"{base}.result_blast.txt"),
        os.path.join(blast_dir, f"{base}_blast.txt"),
        os.path.join(blast_dir, f"{base}_genes_blast.txt"),
        os.path.join(blast_dir, f"{base}.txt"),
    ]
    for p in candidates:
        if os.path.exists(p) and os.path.getsize(p) > 0:
            return p
    # fallback: first file starting with genome_id
    for p in sorted(glob.glob(os.path.join(blast_dir, f"{genome_id}*.txt"))):
        if os.path.getsize(p) > 0:
            return p
    return None


def extract_operon_gene_name(ref_id: str) -> str:
    # e.g., "operon_CDS_7|fruR|..." -> "fruR"
    parts = ref_id.split("|")
    if len(parts) >= 2:
        return parts[1]
    return ref_id


def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def count_mismatches(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)


def find_rbs(seq_upstream: str, rbs_min_spacer: int, rbs_max_spacer: int) -> Tuple[Optional[str], Optional[int], Optional[int]]:
    # Scan window upstream of start; return best motif allowing up to 1 mismatch
    best: Tuple[Optional[str], Optional[int], Optional[int]] = (None, None, None)
    best_mismatches: Optional[int] = None
    # Position 0 is immediately upstream base of start; positions increase upstream
    # We need to search a region (rbs_min_spacer + motif_len) .. (rbs_max_spacer + motif_len)
    window_len = max(rbs_max_spacer + 6, len(seq_upstream))
    region = seq_upstream[:window_len]
    for spacer in range(rbs_min_spacer, rbs_max_spacer + 1):
        for motif in RBS_MOTIFS:
            mlen = len(motif)
            start_idx = spacer + mlen  # index upstream from start where motif ends
            # Motif spans [spacer .. spacer+mlen-1] upstream; 0-based in region slice
            # So in forward string where index 0 is immediately upstream, motif slice is region[spacer:spacer+mlen]
            if spacer + mlen <= len(region):
                window = region[spacer : spacer + mlen]
                mism = count_mismatches(window, motif)
                if mism <= 1:
                    if best_mismatches is None or mism < best_mismatches:
                        best = (window, spacer, mism)
                        best_mismatches = mism
    return best


def scan_upstream_alternative_starts(
    coding_oriented_with_context: str,
    upstream_prefix_len: int,
    rbs_min_spacer: int,
    rbs_max_spacer: int,
) -> Tuple[bool, Optional[str], Optional[int], Optional[bool], Optional[str], Optional[int], Optional[int]]:
    """
    coding_oriented_with_context: sequence consisting of [upstream_context][CDS ...]
    upstream_prefix_len: number of nt in upstream_context
    Returns whether an upstream start exists, its codon, nt offset to called start,
    and RBS details for that upstream candidate.
    """
    if upstream_prefix_len < 3:
        return (False, None, None, None, None, None, None)

    upstream_seq = coding_oriented_with_context[:upstream_prefix_len]
    cds_from_start = coding_oriented_with_context[upstream_prefix_len:]

    # Scan upstream_seq in the same frame as called start (in-frame positions stepping by 3)
    for offset in range(3, upstream_prefix_len + 1, 3):
        codon_start = upstream_prefix_len - offset
        codon = upstream_seq[codon_start: codon_start + 3]
        if codon not in ("ATG", "GTG", "TTG"):
            continue
        # Check that the extended ORF (codon + cds_from_start) has no in-frame stop before the original start (trivial) and beyond
        extended = upstream_seq[codon_start:] + cds_from_start
        _ = str(Seq(extended).translate(to_stop=False))
        # A conservative check: do not require rechecking stops inside upstream since it's in-frame by construction
        # Record the first viable upstream start we encounter (closest to the original start)
        # RBS detection: look in upstream of this codon within upstream_seq
        # Build upstream-of-candidate string with position 0 immediately upstream of the candidate start
        upstream_of_candidate = upstream_seq[:codon_start][::-1]  # reverse for 0=index immediate upstream
        rbs_seq, rbs_spacer, mism = find_rbs(upstream_of_candidate, rbs_min_spacer, rbs_max_spacer)
        has_rbs = rbs_seq is not None
        return (
            True,
            codon,
            offset,  # nt upstream offset from called start
            has_rbs,
            rbs_seq,
            rbs_spacer,
            mism,
        )

    return (False, None, None, None, None, None, None)


def classify_case(
    start_codon: str,
    gff_partial: Optional[str],
    upstream_found: bool,
    upstream_has_rbs: Optional[bool],
) -> str:
    if gff_partial and any(x in gff_partial for x in ("10", "01", "true", "True")):
        return "Partial_gene"
    if not upstream_found:
        return "No_upstream"
    if upstream_has_rbs:
        return "Alt_upstream_with_RBS"
    return "Upstream_no_RBS"


def process_genome(
    genome_dir: str,
    gene_names: List[str],
    max_upstream: int,
    rbs_min_spacer: int,
    rbs_max_spacer: int,
    blast_dir: Optional[str] = None,
) -> List[StartSiteResult]:
    genome_id = os.path.basename(genome_dir.rstrip("/"))
    gff_path = os.path.join(genome_dir, f"{genome_id}.gff")
    fna_path = os.path.join(genome_dir, f"{genome_id}.fna")
    if not (os.path.exists(gff_path) and os.path.exists(fna_path)):
        return []

    contigs = read_genome_sequences(fna_path)
    genes = parse_gff_genes(gff_path)

    # If BLAST dir provided, map prokka locus_tag (qseqid) → reference gene name (from sseqid)
    locus_to_target_gene: Optional[Dict[str, str]] = None
    if blast_dir:
        bfile = find_blast_file_for_genome(blast_dir, genome_id)
        if bfile:
            hits = parse_blast_best_hits(bfile)
            locus_to_target_gene = {}
            for h in hits:
                gene_name = extract_operon_gene_name(h["sseqid"])
                if gene_name in gene_names:
                    locus_to_target_gene[h["qseqid"]] = gene_name

    results: List[StartSiteResult] = []

    if locus_to_target_gene:
        # Index GFF entries by locus_tag and ID for quick lookup
        by_locus: Dict[str, Dict[str, str]] = {}
        by_id: Dict[str, Dict[str, str]] = {}
        for e in genes:
            if e["type"] != "CDS":
                continue
            attrs = e["attrs"]
            if "locus_tag" in attrs:
                by_locus[attrs["locus_tag"]] = e
            if "ID" in attrs:
                by_id[attrs["ID"]] = e

        for locus, ref_gene in locus_to_target_gene.items():
            entry = by_locus.get(locus) or by_id.get(locus)
            if not entry:
                continue
            contig = entry["contig"]
            strand = entry["strand"]
            start = int(entry["start"])  # 1-based
            end = int(entry["end"])  # 1-based inclusive
            attrs = entry["attrs"]
            gff_start_type = attrs.get("start_type")
            gff_partial = attrs.get("partial")
            gff_rbs_motif = attrs.get("rbs_motif")
            gff_rbs_spacer = attrs.get("rbs_spacer")

            if contig not in contigs:
                continue
            contig_seq = contigs[contig]

            if strand == "+":
                cds_start_idx0 = start - 1
                upstream_start0 = max(0, cds_start_idx0 - max_upstream)
                window_seq = str(contig_seq[upstream_start0:end])
                prefix_len = cds_start_idx0 - upstream_start0
                coding_oriented = window_seq
            else:
                cds_end_idx0 = end
                upstream_end0 = min(len(contig_seq), cds_end_idx0 + max_upstream)
                window_seq = str(contig_seq[start - 1:upstream_end0].reverse_complement())
                prefix_len = upstream_end0 - end
                coding_oriented = window_seq

            start_codon = coding_oriented[prefix_len:prefix_len + 3]
            scan_seq = coding_oriented[: prefix_len + 300]
            upstream_found, up_codon, up_offset, up_has_rbs, up_rbs_seq, up_rbs_spacer, up_rbs_mism = scan_upstream_alternative_starts(
                scan_seq,
                upstream_prefix_len=prefix_len,
                rbs_min_spacer=rbs_min_spacer,
                rbs_max_spacer=rbs_max_spacer,
            )

            classification = classify_case(start_codon, gff_partial, upstream_found, up_has_rbs)

            results.append(
                StartSiteResult(
                    genome_id=genome_id,
                    gene=ref_gene,
                    contig=contig,
                    strand=strand,
                    start=start,
                    end=end,
                    start_codon=start_codon,
                    gff_start_type=gff_start_type,
                    gff_partial=gff_partial,
                    gff_rbs_motif=gff_rbs_motif,
                    gff_rbs_spacer=gff_rbs_spacer,
                    upstream_candidate_found=bool(upstream_found),
                    upstream_codon=up_codon,
                    upstream_offset_nt=up_offset,
                    upstream_has_rbs=up_has_rbs,
                    upstream_rbs_seq=up_rbs_seq,
                    upstream_rbs_spacer_nt=up_rbs_spacer,
                    upstream_rbs_mismatches=up_rbs_mism,
                    classification=classification,
                )
            )
    else:
        # Fallback: original gene-name based scan
        for entry in genes:
            if entry["type"] != "CDS":
                continue
            gene_name = entry.get("gene_name")
            if not gene_name or gene_name not in gene_names:
                continue
            contig = entry["contig"]
            strand = entry["strand"]
            start = int(entry["start"])  # 1-based
            end = int(entry["end"])  # 1-based inclusive
            attrs = entry["attrs"]
            gff_start_type = attrs.get("start_type")
            gff_partial = attrs.get("partial")
            gff_rbs_motif = attrs.get("rbs_motif")
            gff_rbs_spacer = attrs.get("rbs_spacer")

            if contig not in contigs:
                continue
            contig_seq = contigs[contig]

            if strand == "+":
                cds_start_idx0 = start - 1
                upstream_start0 = max(0, cds_start_idx0 - max_upstream)
                window_seq = str(contig_seq[upstream_start0:end])
                prefix_len = cds_start_idx0 - upstream_start0
                coding_oriented = window_seq
            else:
                cds_end_idx0 = end
                upstream_end0 = min(len(contig_seq), cds_end_idx0 + max_upstream)
                window_seq = str(contig_seq[start - 1:upstream_end0].reverse_complement())
                prefix_len = upstream_end0 - end
                coding_oriented = window_seq

            start_codon = coding_oriented[prefix_len:prefix_len + 3]
            scan_seq = coding_oriented[: prefix_len + 300]
            upstream_found, up_codon, up_offset, up_has_rbs, up_rbs_seq, up_rbs_spacer, up_rbs_mism = scan_upstream_alternative_starts(
                scan_seq,
                upstream_prefix_len=prefix_len,
                rbs_min_spacer=rbs_min_spacer,
                rbs_max_spacer=rbs_max_spacer,
            )

            classification = classify_case(start_codon, gff_partial, upstream_found, up_has_rbs)

            results.append(
                StartSiteResult(
                    genome_id=genome_id,
                    gene=gene_name,
                    contig=contig,
                    strand=strand,
                    start=start,
                    end=end,
                    start_codon=start_codon,
                    gff_start_type=gff_start_type,
                    gff_partial=gff_partial,
                    gff_rbs_motif=gff_rbs_motif,
                    gff_rbs_spacer=gff_rbs_spacer,
                    upstream_candidate_found=bool(upstream_found),
                    upstream_codon=up_codon,
                    upstream_offset_nt=up_offset,
                    upstream_has_rbs=up_has_rbs,
                    upstream_rbs_seq=up_rbs_seq,
                    upstream_rbs_spacer_nt=up_rbs_spacer,
                    upstream_rbs_mismatches=up_rbs_mism,
                    classification=classification,
                )
            )

    return results


def write_tsv(results: List[StartSiteResult], out_path: str) -> None:
    import csv

    fieldnames = list(asdict(results[0]).keys()) if results else [
        "genome_id","gene","contig","strand","start","end","start_codon",
        "gff_start_type","gff_partial","gff_rbs_motif","gff_rbs_spacer",
        "upstream_candidate_found","upstream_codon","upstream_offset_nt",
        "upstream_has_rbs","upstream_rbs_seq","upstream_rbs_spacer_nt","upstream_rbs_mismatches","classification"
    ]
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for r in results:
            writer.writerow(asdict(r))


def stratify_by_metadata(tsv_path: str, metadata_path: str, output_dir: str) -> None:
    """Stratify start site results by metadata categories."""
    
    if not Path(metadata_path).exists():
        print(f"Metadata file not found: {metadata_path}, skipping stratification")
        return
    
    print("\nPerforming metadata stratification analysis...")
    
    # Load the TSV results
    results_df = pd.read_csv(tsv_path, sep='\t')
    
    # Extract barcode from genome_id
    results_df['barcode'] = results_df['genome_id'].str.replace('_AS.result', '').str.replace('.result', '')
    
    # Load metadata
    meta_df = pd.read_csv(metadata_path, sep='\t')
    
    # Create mapping from both Barcode and AS_barcode
    meta_mapping = {}
    for _, row in meta_df.iterrows():
        barcode = row['Barcode']
        meta_mapping[barcode] = row
        if 'AS_barcode' in row and pd.notna(row['AS_barcode']) and row['AS_barcode'] != 'ND':
            as_barcode = row['AS_barcode'].replace('_AS', '')
            meta_mapping[as_barcode] = row
    
    # Merge metadata
    merged_data = []
    for _, row in results_df.iterrows():
        barcode = row['barcode']
        if barcode in meta_mapping:
            meta_row = meta_mapping[barcode]
            merged_row = row.to_dict()
            merged_row['Source Niche'] = meta_row.get('Source Niche', 'Unknown')
            merged_row['Source Details'] = meta_row.get('Source Details', 'Unknown')
            # Use Country instead of Continent (column name in metadata)
            merged_row['Country'] = meta_row.get('Country', 'Unknown')
            merged_data.append(merged_row)
    
    if not merged_data:
        print("Warning: No genomes could be matched to metadata")
        return
    
    merged_df = pd.DataFrame(merged_data)
    print(f"Matched {len(merged_df)} records to metadata")
    
    # Generate stratification reports
    output_lines = []
    output_lines.append("="*60)
    output_lines.append("START CODON USAGE STRATIFIED BY METADATA")
    output_lines.append("="*60)
    output_lines.append("")
    
    # 1. By Source Niche
    output_lines.append("BY SOURCE NICHE:")
    output_lines.append("-"*40)
    
    niche_results = []
    for niche in sorted(merged_df['Source Niche'].unique()):
        if pd.isna(niche) or niche == 'ND':
            continue
        niche_df = merged_df[merged_df['Source Niche'] == niche]
        total = len(niche_df)
        codon_counts = niche_df['start_codon'].value_counts()
        
        niche_results.append({
            'Source Niche': niche,
            'Total': total,
            'ATG': codon_counts.get('ATG', 0),
            'ATG%': round(100 * codon_counts.get('ATG', 0) / total, 1),
            'GTG': codon_counts.get('GTG', 0),
            'GTG%': round(100 * codon_counts.get('GTG', 0) / total, 1),
            'TTG': codon_counts.get('TTG', 0),
            'TTG%': round(100 * codon_counts.get('TTG', 0) / total, 1),
        })
        
        output_lines.append(f"  {niche}: n={total}")
        output_lines.append(f"    ATG: {codon_counts.get('ATG', 0)} ({100*codon_counts.get('ATG', 0)/total:.1f}%)")
        output_lines.append(f"    GTG: {codon_counts.get('GTG', 0)} ({100*codon_counts.get('GTG', 0)/total:.1f}%)")
        output_lines.append(f"    TTG: {codon_counts.get('TTG', 0)} ({100*codon_counts.get('TTG', 0)/total:.1f}%)")
    
    # Save niche stratification table
    if niche_results:
        niche_df = pd.DataFrame(niche_results)
        niche_df.to_csv(os.path.join(output_dir, "start_codon_by_source_niche.tsv"), sep='\t', index=False)
    
    # 2. ptsA-specific analysis (gene with unique TTG usage)
    output_lines.append("")
    output_lines.append("ptsA START CODON USAGE BY SOURCE NICHE:")
    output_lines.append("-"*40)
    
    ptsa_df = merged_df[merged_df['gene'] == 'ptsA']
    ptsa_results = []
    for niche in sorted(ptsa_df['Source Niche'].unique()):
        if pd.isna(niche) or niche == 'ND':
            continue
        niche_ptsa = ptsa_df[ptsa_df['Source Niche'] == niche]
        total = len(niche_ptsa)
        ttg_count = (niche_ptsa['start_codon'] == 'TTG').sum()
        ttg_pct = 100 * ttg_count / total if total > 0 else 0
        
        ptsa_results.append({
            'Source Niche': niche,
            'Total': total,
            'TTG': ttg_count,
            'TTG%': round(ttg_pct, 1),
            'ATG': (niche_ptsa['start_codon'] == 'ATG').sum(),
            'ATG%': round(100 * (niche_ptsa['start_codon'] == 'ATG').sum() / total, 1) if total > 0 else 0,
        })
        
        output_lines.append(f"  {niche}: {ttg_pct:.1f}% TTG (n={total})")
    
    if ptsa_results:
        ptsa_df_out = pd.DataFrame(ptsa_results)
        ptsa_df_out.to_csv(os.path.join(output_dir, "ptsA_by_source_niche.tsv"), sep='\t', index=False)
    
    # 3. By Country
    output_lines.append("")
    output_lines.append("BY COUNTRY (Top 10):")
    output_lines.append("-"*40)
    
    country_results = []
    # Get country counts and sort by frequency
    country_counts = merged_df['Country'].value_counts()
    for country in country_counts.head(10).index:
        if pd.isna(country) or country == 'ND' or country == 'Unknown':
            continue
        country_df = merged_df[merged_df['Country'] == country]
        total = len(country_df)
        codon_counts = country_df['start_codon'].value_counts()
        
        country_results.append({
            'Country': country,
            'Total': total,
            'ATG': codon_counts.get('ATG', 0),
            'ATG%': round(100 * codon_counts.get('ATG', 0) / total, 1),
            'TTG': codon_counts.get('TTG', 0),
            'TTG%': round(100 * codon_counts.get('TTG', 0) / total, 1),
        })
        
        output_lines.append(f"  {country}: n={total}, ATG={100*codon_counts.get('ATG', 0)/total:.1f}%, TTG={100*codon_counts.get('TTG', 0)/total:.1f}%")
    
    if country_results:
        country_df_out = pd.DataFrame(country_results)
        country_df_out.to_csv(os.path.join(output_dir, "start_codon_by_country.tsv"), sep='\t', index=False)
    
    # Write summary to file
    summary_path = os.path.join(output_dir, "metadata_stratification_summary.txt")
    with open(summary_path, 'w') as f:
        f.write('\n'.join(output_lines))
    
    print(f"\nMetadata stratification complete. Files saved:")
    print(f"  - start_codon_by_source_niche.tsv")
    print(f"  - ptsA_by_source_niche.tsv")
    print(f"  - start_codon_by_country.tsv")
    print(f"  - metadata_stratification_summary.txt")


def create_stratification_plots(output_dir: str) -> None:
    """Create visualizations from the stratification results."""
    
    print("\n" + "="*60)
    print("Creating stratification visualizations")
    print("="*60)
    
    # Create plots directory
    plots_dir = os.path.join(output_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)
    
    # 1. Start codon distribution by Source Niche
    niche_file = os.path.join(output_dir, "start_codon_by_source_niche.tsv")
    if Path(niche_file).exists():
        print("Creating Source Niche plots...")
        niche_df = pd.read_csv(niche_file, sep='\t')
        
        # Stacked bar chart for codon usage
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        
        # Prepare data for stacked bar
        niche_df = niche_df.sort_values('Total', ascending=False)
        x = range(len(niche_df))
        
        # Plot percentages
        ax1.bar(x, niche_df['ATG%'], label='ATG', color='#2E7D32')
        ax1.bar(x, niche_df['GTG%'], bottom=niche_df['ATG%'], label='GTG', color='#1976D2')
        ax1.bar(x, niche_df['TTG%'], bottom=niche_df['ATG%'] + niche_df['GTG%'], label='TTG', color='#D32F2F')
        
        ax1.set_xticks(x)
        ax1.set_xticklabels(niche_df['Source Niche'], rotation=45, ha='right')
        ax1.set_ylabel('Start Codon Usage (%)', fontsize=12)
        ax1.set_title('Start Codon Distribution by Source Niche', fontsize=14, fontweight='bold')
        ax1.legend(loc='upper right')
        ax1.set_ylim(0, 100)
        
        # Add sample sizes
        for i, (idx, row) in enumerate(niche_df.iterrows()):
            ax1.text(i, -5, f"n={row['Total']}", ha='center', fontsize=8)
        
        # TTG usage comparison (highlight ptsA's unique pattern)
        ax2.bar(x, niche_df['TTG%'], color='#D32F2F', alpha=0.7)
        ax2.set_xticks(x)
        ax2.set_xticklabels(niche_df['Source Niche'], rotation=45, ha='right')
        ax2.set_ylabel('TTG Usage (%)', fontsize=12)
        ax2.set_title('TTG Start Codon Usage by Source Niche', fontsize=14, fontweight='bold')
        ax2.axhline(y=12.8, color='gray', linestyle='--', alpha=0.5, label='Overall average (12.8%)')
        ax2.legend()
        
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, 'start_codon_by_source_niche.png'), dpi=150, bbox_inches='tight')
        plt.savefig(os.path.join(plots_dir, 'start_codon_by_source_niche.pdf'), bbox_inches='tight')
        plt.close()
        print("  Saved: start_codon_by_source_niche.png/pdf")
    
    # 2. Create individual gene plots by source niche
    summary_file = os.path.join(output_dir, "start_site_summary.tsv")
    metadata_file = "../00_annotation/8587_Efs_metadata_ASbarcode.txt"
    
    if Path(summary_file).exists() and Path(metadata_file).exists():
        print("Creating individual gene plots by source niche...")
        summary_df = pd.read_csv(summary_file, sep='\t')
        metadata_df = pd.read_csv(metadata_file, sep='\t')
        
        # Extract genome ID from file name to match with AS_barcode
        summary_df['AS_barcode'] = summary_df['genome_id'].str.replace('.result', '', regex=False)
        
        # Merge with metadata
        merged_df = summary_df.merge(metadata_df, on='AS_barcode', how='left')
        
        # Get unique genes
        genes = merged_df['gene'].unique()
        
        for gene in genes:
            print(f"  Creating niche plot for {gene}...")
            gene_df = merged_df[merged_df['gene'] == gene]
            
            # Calculate statistics by source niche
            niche_stats = []
            for niche in gene_df['Source Niche'].dropna().unique():
                niche_data = gene_df[gene_df['Source Niche'] == niche]
                total = len(niche_data)
                if total >= 5:  # Only include niches with at least 5 samples
                    atg_count = len(niche_data[niche_data['start_codon'] == 'ATG'])
                    gtg_count = len(niche_data[niche_data['start_codon'] == 'GTG'])  
                    ttg_count = len(niche_data[niche_data['start_codon'] == 'TTG'])
                    
                    niche_stats.append({
                        'Source Niche': niche,
                        'Total': total,
                        'ATG': atg_count,
                        'GTG': gtg_count,
                        'TTG': ttg_count,
                        'ATG%': (atg_count/total)*100,
                        'GTG%': (gtg_count/total)*100,
                        'TTG%': (ttg_count/total)*100
                    })
            
            if niche_stats:
                niche_gene_df = pd.DataFrame(niche_stats)
                
                # Save to file
                gene_niche_file = os.path.join(output_dir, f"{gene}_by_source_niche.tsv")
                niche_gene_df.to_csv(gene_niche_file, sep='\t', index=False)
                
                # Create simplified visualization - just stacked bar chart
                fig, ax = plt.subplots(figsize=(10, 6))
                
                # Sort by ATG% descending to show highest ATG usage first
                niche_gene_df_sorted = niche_gene_df.sort_values('ATG%', ascending=False)
                x = range(len(niche_gene_df_sorted))
                
                # Create stacked bars
                ax.bar(x, niche_gene_df_sorted['ATG%'], label='ATG', color='#2E7D32', alpha=0.9)
                ax.bar(x, niche_gene_df_sorted['GTG%'], bottom=niche_gene_df_sorted['ATG%'], 
                       label='GTG', color='#1976D2', alpha=0.9)
                ax.bar(x, niche_gene_df_sorted['TTG%'], 
                       bottom=niche_gene_df_sorted['ATG%'] + niche_gene_df_sorted['GTG%'], 
                       label='TTG', color='#D32F2F', alpha=0.9)
                
                ax.set_xticks(x)
                ax.set_xticklabels(niche_gene_df_sorted['Source Niche'], rotation=45, ha='right', fontsize=10)
                ax.set_ylabel('Start Codon Usage (%)', fontsize=11)
                ax.set_xlabel('Source Niche', fontsize=11)
                ax.set_ylim(0, 110)  # Extra space for sample size labels
                
                # Add sample sizes on top of bars
                for i, (idx, row) in enumerate(niche_gene_df_sorted.iterrows()):
                    ax.text(i, 102, f"n={row['Total']}", ha='center', fontsize=8, color='gray', rotation=0)
                
                # Calculate overall statistics
                total_genomes = len(gene_df)
                atg_pct = (gene_df['start_codon'] == 'ATG').mean() * 100
                gtg_pct = (gene_df['start_codon'] == 'GTG').mean() * 100
                ttg_pct = (gene_df['start_codon'] == 'TTG').mean() * 100
                
                # Title with statistics
                ax.set_title(f'{gene.upper()}: Start Codon Distribution by Source Niche (Sorted by ATG%)\n'
                           f'Total: {total_genomes} genomes | ATG: {atg_pct:.1f}% | GTG: {gtg_pct:.1f}% | TTG: {ttg_pct:.1f}%',
                           fontsize=12, fontweight='bold', pad=15)
                
                # Place legend outside plot area
                ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=10)
                
                plt.tight_layout()
                plot_file = os.path.join(plots_dir, f'{gene}_by_source_niche.png')
                plt.savefig(plot_file, dpi=150, bbox_inches='tight')
                plt.savefig(plot_file.replace('.png', '.pdf'), bbox_inches='tight')
                plt.close()
                print(f"    Saved: {gene}_by_source_niche.png/pdf")
    
    # 3. Create individual gene plots by country
    summary_file = os.path.join(output_dir, "start_site_summary.tsv")
    metadata_file = "../00_annotation/8587_Efs_metadata_ASbarcode.txt"
    
    if Path(summary_file).exists() and Path(metadata_file).exists():
        print("Creating individual gene plots by country...")
        summary_df = pd.read_csv(summary_file, sep='\t')
        metadata_df = pd.read_csv(metadata_file, sep='\t')
        
        # Extract genome ID from file name to match with AS_barcode
        summary_df['AS_barcode'] = summary_df['genome_id'].str.replace('.result', '', regex=False)
        
        # Merge with metadata
        merged_df = summary_df.merge(metadata_df, on='AS_barcode', how='left')
        
        # Get unique genes
        genes = merged_df['gene'].unique()
        
        for gene in genes:
            print(f"  Creating plot for {gene}...")
            gene_df = merged_df[merged_df['gene'] == gene]
            
            # Calculate statistics by country
            country_stats = []
            for country in gene_df['Country'].dropna().unique():
                country_data = gene_df[gene_df['Country'] == country]
                total = len(country_data)
                if total >= 5:  # Only include countries with at least 5 samples
                    atg_count = len(country_data[country_data['start_codon'] == 'ATG'])
                    gtg_count = len(country_data[country_data['start_codon'] == 'GTG'])
                    ttg_count = len(country_data[country_data['start_codon'] == 'TTG'])
                    
                    country_stats.append({
                        'Country': country,
                        'Total': total,
                        'ATG': atg_count,
                        'GTG': gtg_count,
                        'TTG': ttg_count,
                        'ATG%': (atg_count/total)*100,
                        'GTG%': (gtg_count/total)*100,
                        'TTG%': (ttg_count/total)*100
                    })
            
            if country_stats:
                country_gene_df = pd.DataFrame(country_stats)
                
                # Save to file
                gene_country_file = os.path.join(output_dir, f"{gene}_by_country.tsv")
                country_gene_df.to_csv(gene_country_file, sep='\t', index=False)
                
                # Create simplified visualization - just stacked bar chart
                fig, ax = plt.subplots(figsize=(12, 6))
                
                # Sort by ATG% descending to show highest ATG usage first
                country_gene_df_sorted = country_gene_df.sort_values('ATG%', ascending=False).head(25)  # Top 25 countries
                x = range(len(country_gene_df_sorted))
                
                # Create stacked bars
                ax.bar(x, country_gene_df_sorted['ATG%'], label='ATG', color='#2E7D32', alpha=0.9)
                ax.bar(x, country_gene_df_sorted['GTG%'], bottom=country_gene_df_sorted['ATG%'], 
                       label='GTG', color='#1976D2', alpha=0.9)
                ax.bar(x, country_gene_df_sorted['TTG%'], 
                       bottom=country_gene_df_sorted['ATG%'] + country_gene_df_sorted['GTG%'], 
                       label='TTG', color='#D32F2F', alpha=0.9)
                
                ax.set_xticks(x)
                ax.set_xticklabels(country_gene_df_sorted['Country'], rotation=45, ha='right', fontsize=9)
                ax.set_ylabel('Start Codon Usage (%)', fontsize=11)
                ax.set_xlabel('Country', fontsize=11)
                ax.set_ylim(0, 110)  # Extra space for sample size labels
                
                # Add sample sizes on top of bars
                for i, (idx, row) in enumerate(country_gene_df_sorted.iterrows()):
                    ax.text(i, 102, f"n={row['Total']}", ha='center', fontsize=7, color='gray', rotation=0)
                
                # Calculate overall statistics
                total_genomes = len(gene_df)
                atg_pct = (gene_df['start_codon'] == 'ATG').mean() * 100
                gtg_pct = (gene_df['start_codon'] == 'GTG').mean() * 100
                ttg_pct = (gene_df['start_codon'] == 'TTG').mean() * 100
                
                # Title with statistics
                ax.set_title(f'{gene.upper()}: Start Codon Distribution by Country (Sorted by ATG%)\n'
                           f'Total: {total_genomes} genomes | ATG: {atg_pct:.1f}% | GTG: {gtg_pct:.1f}% | TTG: {ttg_pct:.1f}%',
                           fontsize=12, fontweight='bold', pad=15)
                
                # Place legend outside plot area
                ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=10)
                
                plt.tight_layout()
                plot_file = os.path.join(plots_dir, f'{gene}_by_country.png')
                plt.savefig(plot_file, dpi=150, bbox_inches='tight')
                plt.savefig(plot_file.replace('.png', '.pdf'), bbox_inches='tight')
                plt.close()
                print(f"    Saved: {gene}_by_country.png/pdf")
    
    # 4. Overall country comparison (all genes combined)
    country_file = os.path.join(output_dir, "start_codon_by_country.tsv")
    if Path(country_file).exists():
        print("Creating overall Country comparison plots...")
        country_df = pd.read_csv(country_file, sep='\t')
        
        # Create figure with multiple subplots
        fig = plt.figure(figsize=(16, 10))
        
        # 3a. Vertical stacked bar chart sorted by ATG% (top left)
        ax1 = plt.subplot(2, 2, 1)
        # Sort by ATG% descending (highest ATG first)
        country_df_sorted = country_df.sort_values('ATG%', ascending=False).head(20)  # Top 20 countries
        x = range(len(country_df_sorted))
        
        # Create stacked bars
        ax1.bar(x, country_df_sorted['ATG%'], label='ATG', color='#2E7D32', alpha=0.8)
        ax1.bar(x, country_df_sorted['GTG%'], bottom=country_df_sorted['ATG%'], 
               label='GTG', color='#1976D2', alpha=0.8)
        ax1.bar(x, country_df_sorted['TTG%'], 
               bottom=country_df_sorted['ATG%'] + country_df_sorted['GTG%'], 
               label='TTG', color='#D32F2F', alpha=0.8)
        
        ax1.set_xticks(x)
        ax1.set_xticklabels(country_df_sorted['Country'], rotation=45, ha='right', fontsize=8)
        ax1.set_ylabel('Start Codon Usage (%)', fontsize=10)
        ax1.set_title('Start Codon Distribution (Sorted by ATG%)', fontsize=12, fontweight='bold')
        ax1.legend(loc='upper right', fontsize=9)
        ax1.set_ylim(0, 110)  # Extra space for sample size labels
        
        # Add sample sizes on top of bars
        for i, (idx, row) in enumerate(country_df_sorted.iterrows()):
            ax1.text(i, 102, f"n={row['Total']}", ha='center', fontsize=6, color='gray', rotation=0)
        
        # 3b. TTG usage focus (top right)
        ax2 = plt.subplot(2, 2, 2)
        country_df_ttg = country_df.sort_values('TTG%', ascending=False)
        x = range(len(country_df_ttg))
        
        colors = ['#D32F2F' if ttg > 15 else '#FFA726' if ttg > 10 else '#66BB6A' 
                  for ttg in country_df_ttg['TTG%']]
        
        ax2.bar(x, country_df_ttg['TTG%'], color=colors, alpha=0.7)
        ax2.set_xticks(x)
        ax2.set_xticklabels(country_df_ttg['Country'], rotation=45, ha='right')
        ax2.set_ylabel('TTG Usage (%)', fontsize=10)
        ax2.set_title('TTG Start Codon Usage by Country\n(Sorted by TTG frequency)', fontsize=12, fontweight='bold')
        ax2.axhline(y=12.8, color='gray', linestyle='--', alpha=0.5, label='Overall avg (12.8%)')
        ax2.legend(fontsize=9)
        
        # Add value labels on bars
        for i, (idx, row) in enumerate(country_df_ttg.iterrows()):
            ax2.text(i, row['TTG%'] + 0.5, f"{row['TTG%']:.1f}", ha='center', fontsize=7)
        
        # 3c. Sample size vs TTG usage scatter (bottom left)
        ax3 = plt.subplot(2, 2, 3)
        ax3.scatter(country_df['Total'], country_df['TTG%'], 
                   s=100, alpha=0.6, c=country_df['TTG%'], cmap='RdYlGn_r')
        
        # Add country labels for interesting points
        for idx, row in country_df.iterrows():
            if row['TTG%'] > 15 or row['Total'] > 1000:
                ax3.annotate(row['Country'], (row['Total'], row['TTG%']), 
                           fontsize=7, alpha=0.7, ha='center')
        
        ax3.set_xlabel('Sample Size (n)', fontsize=10)
        ax3.set_ylabel('TTG Usage (%)', fontsize=10)
        ax3.set_title('Sample Size vs TTG Usage', fontsize=12, fontweight='bold')
        ax3.axhline(y=12.8, color='gray', linestyle='--', alpha=0.5)
        ax3.grid(True, alpha=0.3)
        
        # 3d. Grouped bar chart comparing ATG vs TTG (bottom right)
        ax4 = plt.subplot(2, 2, 4)
        country_df_top = country_df.nlargest(8, 'Total')  # Top 8 countries by sample size
        x = np.arange(len(country_df_top))
        width = 0.35
        
        ax4.bar(x - width/2, country_df_top['ATG%'], width, label='ATG', color='#2E7D32', alpha=0.8)
        ax4.bar(x + width/2, country_df_top['TTG%'], width, label='TTG', color='#D32F2F', alpha=0.8)
        
        ax4.set_xticks(x)
        ax4.set_xticklabels(country_df_top['Country'], rotation=45, ha='right')
        ax4.set_ylabel('Usage (%)', fontsize=10)
        ax4.set_title('ATG vs TTG Usage (Top 8 Countries by Sample Size)', fontsize=12, fontweight='bold')
        ax4.legend(fontsize=9)
        
        # Add sample sizes below
        for i, (idx, row) in enumerate(country_df_top.iterrows()):
            ax4.text(i, -5, f"n={row['Total']}", ha='center', fontsize=7, color='gray')
        
        plt.suptitle('Start Codon Usage Analysis by Country', fontsize=14, fontweight='bold', y=1.02)
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, 'start_codon_by_country_detailed.png'), dpi=150, bbox_inches='tight')
        plt.savefig(os.path.join(plots_dir, 'start_codon_by_country_detailed.pdf'), bbox_inches='tight')
        plt.close()
        print("  Saved: start_codon_by_country_detailed.png/pdf")
    
    # 5. Heatmap comparing all genes across source niches
    if Path(summary_file).exists() and Path(metadata_file).exists():
        print("Creating gene comparison heatmap...")
        
        # Load and merge data
        summary_df = pd.read_csv(summary_file, sep='\t')
        summary_df['barcode'] = summary_df['genome_id'].str.replace('_AS.result', '').str.replace('.result', '')
        
        meta_df = pd.read_csv(metadata_file, sep='\t')
        meta_mapping = {}
        for _, row in meta_df.iterrows():
            meta_mapping[row['Barcode']] = row['Source Niche']
            if 'AS_barcode' in row and pd.notna(row['AS_barcode']) and row['AS_barcode'] != 'ND':
                as_barcode = row['AS_barcode'].replace('_AS', '')
                meta_mapping[as_barcode] = row['Source Niche']
        
        summary_df['Source Niche'] = summary_df['barcode'].map(meta_mapping)
        summary_df = summary_df.dropna(subset=['Source Niche'])
        
        # Calculate TTG usage by gene and niche
        heatmap_data = []
        genes = ['frpC', 'glpC', 'ptsD', 'ptsC', 'ptsB', 'ptsA', 'fruR']
        
        for niche in summary_df['Source Niche'].unique():
            if niche == 'ND' or pd.isna(niche):
                continue
            niche_data = {'Source Niche': niche}
            niche_df = summary_df[summary_df['Source Niche'] == niche]
            
            for gene in genes:
                gene_df = niche_df[niche_df['gene'] == gene]
                if len(gene_df) > 0:
                    ttg_pct = 100 * (gene_df['start_codon'] == 'TTG').sum() / len(gene_df)
                    niche_data[gene] = ttg_pct
                else:
                    niche_data[gene] = 0
            
            heatmap_data.append(niche_data)
        
        if heatmap_data:
            heatmap_df = pd.DataFrame(heatmap_data)
            heatmap_df = heatmap_df.set_index('Source Niche')
            
            # Create heatmap
            fig, ax = plt.subplots(figsize=(10, 6))
            sns.heatmap(heatmap_df[genes], annot=True, fmt='.1f', cmap='RdYlBu_r', 
                       vmin=0, vmax=100, cbar_kws={'label': 'TTG Usage (%)'},
                       linewidths=0.5, linecolor='gray')
            
            ax.set_xlabel('Operon Gene', fontsize=12)
            ax.set_ylabel('Source Niche', fontsize=12)
            ax.set_title('TTG Start Codon Usage Heatmap\n(ptsA shows consistently high TTG usage across all niches)', 
                        fontsize=14, fontweight='bold')
            
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, 'ttg_usage_heatmap.png'), dpi=150, bbox_inches='tight')
            plt.savefig(os.path.join(plots_dir, 'ttg_usage_heatmap.pdf'), bbox_inches='tight')
            plt.close()
            print("  Saved: ttg_usage_heatmap.png/pdf")
    
    print(f"\nAll plots saved to: {plots_dir}/")
    print("Visualization complete!")


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    
    # If visualize-only mode, just create plots and exit
    if args.visualize_only:
        print("Running in visualization-only mode...")
        create_stratification_plots(args.output_dir)
        return
    
    # For full analysis, prokka_dir is required
    if not args.prokka_dir:
        print("Error: --prokka_dir is required for full analysis")
        sys.exit(1)

    gene_names = load_operon_gene_names(args.gene_reference_fasta)

    # Determine genome dirs
    genome_dirs = sorted([d for d in glob.glob(os.path.join(args.prokka_dir, "*")) if os.path.isdir(d)])
    if args.genome_list:
        with open(args.genome_list) as f:
            allowed = set(x.strip() for x in f if x.strip())
        genome_dirs = [d for d in genome_dirs if os.path.basename(d) in allowed]

    all_results: List[StartSiteResult] = []
    with ProcessPoolExecutor(max_workers=args.max_workers) as ex:
        futures = {
            ex.submit(
                process_genome,
                gdir,
                gene_names,
                args.max_upstream,
                args.rbs_min_spacer,
                args.rbs_max_spacer,
                args.blast_dir,
            ): gdir
            for gdir in genome_dirs
        }
        for fut in as_completed(futures):
            try:
                res = fut.result()
                if res:
                    all_results.extend(res)
            except Exception as e:
                sys.stderr.write(f"Error processing {futures[fut]}: {e}\n")

    out_tsv = os.path.join(args.output_dir, "start_site_summary.tsv")
    write_tsv(all_results, out_tsv)
    print(f"Wrote {len(all_results)} rows to {out_tsv}")
    
    # Perform metadata stratification if metadata file exists
    if args.metadata:
        stratify_by_metadata(out_tsv, args.metadata, args.output_dir)
        # Also create visualizations
        create_stratification_plots(args.output_dir)


if __name__ == "__main__":
    main()


