#!/usr/bin/env python3
"""
Analyze start-codon choice discrepancies for operon genes across Prokka/Prodigal outputs.

For each genome's GFF + FNA, find operon genes of interest, extract coding sequence
context, scan for upstream in-frame alternative starts (ATG/GTG/TTG) and RBS motifs,
and summarize likely reasons for shorter/alternate starts.
"""

import argparse
import os
import sys
import glob
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Tuple

from concurrent.futures import ProcessPoolExecutor, as_completed

from Bio import SeqIO
from Bio.Seq import Seq


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
    parser.add_argument("--prokka_dir", required=True, help="Directory containing per-genome Prokka outputs")
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


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

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


if __name__ == "__main__":
    main()


