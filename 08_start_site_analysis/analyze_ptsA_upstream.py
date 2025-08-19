#!/usr/bin/env python3
"""
ptsA-focused upstream/start analysis:
- Extract standardized upstream windows around ptsA start across genomes
- Align with MAFFT
- Scan RBS motifs (AGGAGG family) and compute spacer distributions
- Produce summary plots and an HTML report
"""

import os
import sys
import glob
import argparse
from typing import Dict, List, Optional, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RBS_MOTIFS = ["AGGAGG", "GGAGG", "AGGA", "GGAG", "GAGG"]


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def parse_blast_best_hits(blast_path: str) -> List[Dict[str, str]]:
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


def extract_operon_gene_name(ref_id: str) -> str:
    parts = ref_id.split("|")
    if len(parts) >= 2:
        return parts[1]
    return ref_id


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
    for p in sorted(glob.glob(os.path.join(blast_dir, f"{base}*.txt"))):
        if os.path.getsize(p) > 0:
            return p
    return None


def read_genome_sequences(fna_path: str) -> Dict[str, Seq]:
    seqs: Dict[str, Seq] = {}
    for rec in SeqIO.parse(fna_path, "fasta"):
        seqs[rec.id] = rec.seq
    return seqs


def parse_gff_cds(gff_path: str) -> List[Dict[str, str]]:
    cds_list: List[Dict[str, str]] = []
    with open(gff_path, "r") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = parts
            if ftype != "CDS":
                continue
            attr_map: Dict[str, str] = {}
            for a in attrs.split(";"):
                if not a:
                    continue
                if "=" in a:
                    k, v = a.split("=", 1)
                    attr_map[k] = v
            cds_list.append({
                "contig": seqid,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "attrs": attr_map,
            })
    return cds_list


def find_rbs(seq_upstream: str, rbs_min_spacer: int, rbs_max_spacer: int) -> Tuple[Optional[str], Optional[int], Optional[int]]:
    best = (None, None, None)
    best_mism = None
    region = seq_upstream  # string where index 0 is immediately upstream
    for spacer in range(rbs_min_spacer, rbs_max_spacer + 1):
        for motif in RBS_MOTIFS:
            mlen = len(motif)
            if spacer + mlen <= len(region):
                window = region[spacer: spacer + mlen]
                mism = sum(1 for x, y in zip(window, motif) if x != y)
                if mism <= 1 and (best_mism is None or mism < best_mism):
                    best = (window, spacer, mism)
                    best_mism = mism
    return best


def extract_ptsA_windows(
    prokka_dir: str,
    blast_dir: str,
    upstream_len: int,
    downstream_len: int,
    out_fa: str,
    full_gene_downstream: bool = False,
) -> pd.DataFrame:
    rows = []
    records: List[SeqRecord] = []

    genome_dirs = sorted([d for d in glob.glob(os.path.join(prokka_dir, "*")) if os.path.isdir(d)])
    for gdir in genome_dirs:
        genome_id = os.path.basename(gdir.rstrip("/"))
        gff = os.path.join(gdir, f"{genome_id}.gff")
        fna = os.path.join(gdir, f"{genome_id}.fna")
        if not (os.path.exists(gff) and os.path.exists(fna)):
            continue
        bfile = find_blast_file_for_genome(blast_dir, genome_id)
        if not bfile:
            continue
        hits = parse_blast_best_hits(bfile)
        locus_for_ptsA: Optional[str] = None
        for h in hits:
            if extract_operon_gene_name(h["sseqid"]) == "ptsA":
                locus_for_ptsA = h["qseqid"]
                break
        if not locus_for_ptsA:
            continue

        contigs = read_genome_sequences(fna)
        cds_list = parse_gff_cds(gff)
        entry = None
        for e in cds_list:
            attrs = e["attrs"]
            if attrs.get("locus_tag") == locus_for_ptsA or attrs.get("ID") == locus_for_ptsA:
                entry = e
                break
        if not entry:
            continue

        contig = entry["contig"]
        if contig not in contigs:
            continue
        strand = entry["strand"]
        start, end = int(entry["start"]), int(entry["end"])  # 1-based inclusive
        contig_seq = contigs[contig]

        if strand == "+":
            cds_start0 = start - 1
            up0 = max(0, cds_start0 - upstream_len)
            if full_gene_downstream:
                dn_end = min(len(contig_seq), end)  # include full CDS downstream of start
            else:
                dn_end = min(len(contig_seq), start - 1 + downstream_len + 3)
            window = str(contig_seq[up0:dn_end])
            prefix = cds_start0 - up0
        else:
            cds_end0 = end  # for minus strand, called start is at genomic 'end'
            up_end0 = min(len(contig_seq), cds_end0 + upstream_len)
            if full_gene_downstream:
                dn_start0 = max(0, start - 1)  # include full CDS in coding orientation after start
            else:
                dn_start0 = max(0, end - downstream_len - 3)
            window = str(contig_seq[dn_start0:up_end0].reverse_complement())
            prefix = up_end0 - end

        # RBS scan in unaligned sequences
        upstream_of_start = window[:prefix][::-1]
        rbs_seq, rbs_spacer, rbs_mism = find_rbs(upstream_of_start, 4, 13)

        rec_id = f"{genome_id}|{locus_for_ptsA}|{contig}|{strand}|start:{start}|end:{end}|rbs:{rbs_seq or 'NA'}|spacer:{rbs_spacer if rbs_spacer is not None else 'NA'}"
        records.append(SeqRecord(Seq(window), id=rec_id, description=""))
        rows.append({
            "genome_id": genome_id,
            "locus_tag": locus_for_ptsA,
            "contig": contig,
            "strand": strand,
            "start": start,
            "end": end,
            "rbs_seq": rbs_seq,
            "rbs_spacer": rbs_spacer,
            "rbs_mismatches": rbs_mism,
        })

    if records:
        with open(out_fa, "w") as f:
            SeqIO.write(records, f, "fasta")

    return pd.DataFrame(rows)


def run_mafft(input_fa: str, output_fa: str, threads: int = 4) -> bool:
    try:
        # Require MAFFT_TMPDIR to be set and writable
        mafft_tmpdir = os.environ.get("MAFFT_TMPDIR")
        if not mafft_tmpdir:
            raise RuntimeError("MAFFT_TMPDIR is not set")
        if not (os.path.isdir(mafft_tmpdir) and os.access(mafft_tmpdir, os.W_OK)):
            raise RuntimeError(f"MAFFT_TMPDIR is not a writable directory: {mafft_tmpdir}")
        cmd = ["mafft", "--auto", "--thread", str(threads), input_fa]
        with open(output_fa, "w") as out:
            subprocess.run(cmd, check=True, stdout=out, stderr=subprocess.DEVNULL)
        return True
    except Exception as e:
        sys.stderr.write(f"MAFFT failed: {e}\n")
        return False


def plot_rbs_spacer_hist(df: pd.DataFrame, out_png: str) -> None:
    d = df[df["rbs_spacer"].notna()]
    if d.empty:
        return
    plt.figure(figsize=(8,4))
    plt.hist(d["rbs_spacer"], bins=range(3, 15), color="#4C72B0", edgecolor="white")
    plt.xlabel("RBS spacer (nt)")
    plt.ylabel("Count")
    plt.title("ptsA RBS spacer distribution")
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()


def plot_start_codon_counts(summary_tsv: str, out_png: str) -> None:
    df = pd.read_csv(summary_tsv, sep="\t")
    df = df[df["gene"] == "ptsA"]
    if df.empty:
        return
    counts = df["start_codon"].value_counts().sort_values(ascending=False)
    plt.figure(figsize=(6,4))
    counts.plot(kind="bar", color="#55A868")
    plt.ylabel("Genomes")
    plt.title("ptsA called start codons")
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()


def write_html(report_dir: str, assets: Dict[str, str]) -> str:
    html = os.path.join(report_dir, "ptsA_report.html")
    with open(html, "w") as f:
        f.write("<html><head><meta charset='utf-8'><title>ptsA upstream analysis</title></head><body>")
        f.write("<h2>ptsA upstream/start analysis</h2>")
        if os.path.exists(assets.get("start_codon_counts", "")):
            f.write(f"<h3>Called start codons</h3><img src='{os.path.basename(assets['start_codon_counts'])}'><br/>")
        if os.path.exists(assets.get("rbs_spacer_hist", "")):
            f.write(f"<h3>RBS spacer distribution</h3><img src='{os.path.basename(assets['rbs_spacer_hist'])}'><br/>")
        if os.path.exists(assets.get("aligned_fa", "")):
            f.write(f"<h3>Alignment</h3><p>{os.path.basename(assets['aligned_fa'])}</p>")
        if os.path.exists(assets.get("raw_windows", "")):
            f.write(f"<h3>Windows FASTA</h3><p>{os.path.basename(assets['raw_windows'])}</p>")
        f.write("</body></html>")
    return html


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prokka_dir", default="../01_prokka_annotation/output/prokka_results")
    ap.add_argument("--blast_dir", default="../03_blast_search/output/blast_results_prokka_variants")
    ap.add_argument("--summary_tsv", default="output/start_site_summary.tsv")
    ap.add_argument("--out_dir", default="output/ptsA_upstream")
    ap.add_argument("--upstream_len", type=int, default=120)
    ap.add_argument("--downstream_len", type=int, default=60)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--full_gene_downstream", action="store_true", help="If set, extract upstream_len upstream plus the entire CDS downstream of the start site (instead of a fixed downstream_len).")
    args = ap.parse_args()

    ensure_dir(args.out_dir)

    # Extract windows
    raw_fa = os.path.join(args.out_dir, "ptsA_windows.fa")
    meta = extract_ptsA_windows(
        args.prokka_dir,
        args.blast_dir,
        args.upstream_len,
        args.downstream_len,
        raw_fa,
        full_gene_downstream=args.full_gene_downstream,
    )
    meta_path = os.path.join(args.out_dir, "ptsA_windows.tsv")
    meta.to_csv(meta_path, sep="\t", index=False)

    # Align
    aln_fa = os.path.join(args.out_dir, "ptsA_windows.aln.fa")
    aligned = False
    if os.path.exists(raw_fa) and os.path.getsize(raw_fa) > 0:
        aligned = run_mafft(raw_fa, aln_fa, threads=args.threads)

    # Plots
    assets: Dict[str, str] = {}
    sc_png = os.path.join(args.out_dir, "ptsA_start_codon_counts.png")
    plot_start_codon_counts(args.summary_tsv, sc_png)
    assets["start_codon_counts"] = sc_png

    spacer_png = os.path.join(args.out_dir, "ptsA_rbs_spacer_hist.png")
    plot_rbs_spacer_hist(meta, spacer_png)
    assets["rbs_spacer_hist"] = spacer_png

    assets["aligned_fa"] = aln_fa if aligned else ""
    assets["raw_windows"] = raw_fa

    html = write_html(args.out_dir, assets)
    print(f"Wrote ptsA report: {html}")


if __name__ == "__main__":
    main()


