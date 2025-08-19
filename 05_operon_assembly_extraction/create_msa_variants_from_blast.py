#!/usr/bin/env python3
"""
Create MSA(s) using BLAST-derived variant fragments (qseq) to reproduce old alignment behavior.

Inputs:
- blast_dir: directory with per-genome BLAST outputs (qseq in outfmt)
  (e.g., produced by run_blast_nt_vs_prokka_genes.sh)

Outputs:
- output_dir/gene_variants/*.fasta           (per-gene fragments)
- output_dir/msa_variants/*.fasta            (MAFFT alignments of fragments)
"""

import os
import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import os


def parse_blast_qseq_lines(blast_file):
    per_gene = defaultdict(list)
    with open(blast_file, 'r') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            # Expect qseq and sseq present (old-style outfmt)
            if len(parts) < 14:
                continue
            qseqid, sseqid = parts[0], parts[1]
            pident = float(parts[2])
            qseq = parts[12]
            # Derive gene name from subject id (prefer middle token like 'ptsA' if present)
            if '|' in sseqid:
                tokens = sseqid.split('|')
                gene_name = tokens[1] if len(tokens) >= 2 else tokens[-1]
            else:
                gene_name = sseqid
            per_gene[gene_name].append((qseqid, pident, qseq))
    return per_gene


def write_gene_variant_fastas(per_gene, out_dir, min_identity=0.0):
    os.makedirs(out_dir, exist_ok=True)
    written = {}
    for gene, entries in per_gene.items():
        records = []
        for qseqid, pident, qseq in entries:
            if pident < min_identity:
                continue
            seq = qseq.replace('-', '')
            rec = SeqRecord(Seq(seq), id=f"{qseqid}", description=f"pident={pident:.1f}")
            records.append(rec)
        if records:
            out_fa = os.path.join(out_dir, f"{gene}_variants.fasta")
            with open(out_fa, 'w') as fh:
                SeqIO.write(records, fh, 'fasta')
            written[gene] = out_fa
    return written


def run_mafft(input_fa, output_fa, threads=1):
    cmd = ['mafft', '--auto', '--thread', str(threads), input_fa]
    with open(output_fa, 'w') as out:
        env = os.environ.copy()
        env.setdefault('MAFFT_TMPDIR', '/vol/tmp')
        env.setdefault('TMPDIR', env['MAFFT_TMPDIR'])
        res = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True, env=env)
    return res.returncode == 0


def align_variants(variant_fastas, out_dir, threads=1):
    os.makedirs(out_dir, exist_ok=True)
    outputs = []
    for gene, fa in sorted(variant_fastas.items()):
        out_fa = os.path.join(out_dir, f"{gene}_variants_aligned.fasta")
        if run_mafft(fa, out_fa, threads):
            outputs.append(out_fa)
    return outputs


def main():
    ap = argparse.ArgumentParser(description='Create variant-based MSAs from BLAST qseq (old behavior).')
    ap.add_argument('--blast-dir', default='output/blast_nt_vs_prokka', help='Directory with per-genome BLAST txt files')
    ap.add_argument('--output-dir', default='output/variants_pipeline', help='Output base directory')
    ap.add_argument('--min-identity', type=float, default=0.0, help='Minimum %identity to include a variant')
    ap.add_argument('--threads', type=int, default=4)
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    variants_dir = os.path.join(args.output_dir, 'gene_variants')
    msa_dir = os.path.join(args.output_dir, 'msa_variants')

    # Aggregate all per-genome BLAST files
    blast_files = [
        os.path.join(args.blast_dir, f) for f in os.listdir(args.blast_dir)
        if f.endswith('.txt') and os.path.getsize(os.path.join(args.blast_dir, f)) > 0
    ]

    per_gene_all = defaultdict(list)
    for bf in blast_files:
        per_gene = parse_blast_qseq_lines(bf)
        for gene, entries in per_gene.items():
            per_gene_all[gene].extend(entries)

    written = write_gene_variant_fastas(per_gene_all, variants_dir, args.min_identity)
    alignments = align_variants(written, msa_dir, args.threads)

    print(f"Wrote {len(written)} gene variant FASTAs → {variants_dir}")
    print(f"Created {len(alignments)} variant MSAs → {msa_dir}")


if __name__ == '__main__':
    main()


