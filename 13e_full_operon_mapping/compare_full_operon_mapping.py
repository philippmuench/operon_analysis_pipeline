#!/usr/bin/env python3
"""Compare full-operon BLAST mapping for legacy vs SNP-aware GenBank sequences.

For each provided assembly the script runs `blastn` with the full ORIGIN sequence
extracted from two GenBank files (legacy/new) and records the best hit metrics.
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd
from Bio import SeqIO


@dataclass
class BlastHit:
    sseqid: str
    pident: float
    length: int
    evalue: float
    bitscore: float
    qstart: int
    qend: int
    sstart: int
    send: int
    qlen: int
    qseq: str
    sseq: str

    @property
    def coverage(self) -> float:
        if self.qlen == 0:
            return 0.0
        return (self.length / self.qlen) * 100.0


@dataclass
class MappingResult:
    assembly: Path
    legacy_hit: Optional[BlastHit]
    updated_hit: Optional[BlastHit]


def alignment_stats(hit: Optional[BlastHit]) -> Dict[str, Optional[str]]:
    if not hit:
        return {
            "mismatch_count": None,
            "insertion_count": None,
            "deletion_count": None,
            "mismatch_positions": None,
        }

    mismatches: List[str] = []
    mismatch_count = 0
    insertion_count = 0
    deletion_count = 0

    query_pos = hit.qstart
    step = 1 if hit.qend >= hit.qstart else -1

    for q_char, s_char in zip(hit.qseq.upper(), hit.sseq.upper()):
        current_pos = None
        if q_char != "-":
            current_pos = query_pos
            query_pos += step

        if q_char == "-" and s_char != "-":
            insertion_count += 1
        elif q_char != "-" and s_char == "-":
            deletion_count += 1
        elif q_char != "-" and s_char != "-" and q_char != s_char:
            mismatch_count += 1
            if current_pos is not None:
                mismatches.append(f"{current_pos}:{q_char}>{s_char}")

    mismatch_str = ",".join(mismatches[:50])
    if mismatch_count > 50:
        mismatch_str += ",..."

    return {
        "mismatch_count": mismatch_count,
        "insertion_count": insertion_count,
        "deletion_count": deletion_count,
        "mismatch_positions": mismatch_str,
    }


def extract_origin_sequence(gb_path: Path) -> str:
    record = next(SeqIO.parse(str(gb_path), "genbank"))
    return str(record.seq)


def write_temp_fasta(sequence: str, label: str) -> Path:
    tmp = tempfile.NamedTemporaryFile("w", delete=False, suffix=".fasta")
    tmp_path = Path(tmp.name)
    tmp.write(f">{label}\n")
    for i in range(0, len(sequence), 80):
        tmp.write(sequence[i : i + 80] + "\n")
    tmp.close()
    return tmp_path


def ensure_uncompressed(path: Path) -> Path:
    if path.suffix == ".gz":
        tmp = tempfile.NamedTemporaryFile("wb", delete=False, suffix=".fasta")
        with gzip.open(path, "rb") as src, open(tmp.name, "wb") as dst:
            shutil.copyfileobj(src, dst)
        return Path(tmp.name)
    return path


def run_blast(query: Path, subject: Path, threads: int) -> List[BlastHit]:
    subject_path = ensure_uncompressed(subject)
    cmd = [
        "blastn",
        "-query",
        str(query),
        "-subject",
        str(subject_path),
        "-task",
        "megablast",
        "-outfmt",
        "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qseq sseq",
        "-num_threads",
        str(threads),
    ]
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    finally:
        if subject_path is not subject:
            subject_path.unlink(missing_ok=True)
    hits: List[BlastHit] = []
    for line in result.stdout.strip().splitlines():
        parts = line.split("\t")
        if len(parts) < 15:
            continue
        hits.append(
            BlastHit(
                sseqid=parts[1],
                pident=float(parts[2]),
                length=int(parts[3]),
                evalue=float(parts[10]),
                bitscore=float(parts[11]),
                qstart=int(parts[6]),
                qend=int(parts[7]),
                sstart=int(parts[8]),
                send=int(parts[9]),
                qlen=int(parts[12]),
                qseq=parts[13],
                sseq=parts[14],
            )
        )
    hits.sort(key=lambda h: h.bitscore, reverse=True)
    return hits


def summarise_results(results: List[MappingResult]) -> Dict[str, float]:
    summary = {"total": len(results), "legacy_hits": 0, "updated_hits": 0, "legacy_hq": 0, "updated_hq": 0}
    for result in results:
        if result.legacy_hit:
            summary["legacy_hits"] += 1
            if result.legacy_hit.coverage >= 95 and result.legacy_hit.pident >= 95:
                summary["legacy_hq"] += 1
        if result.updated_hit:
            summary["updated_hits"] += 1
            if result.updated_hit.coverage >= 95 and result.updated_hit.pident >= 95:
                summary["updated_hq"] += 1
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-gb", type=Path, default=Path("../02_reference_operon_extraction/operon.gb"))
    parser.add_argument("--updated-gb", type=Path, default=Path("../13a_new_reference_operon_extraction/FL_operon_with_SNPs.gb"))
    parser.add_argument("--assemblies-dir", type=Path, default=Path("../../Efs_assemblies"))
    parser.add_argument("--output", type=Path, default=Path("output"))
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--max-assemblies", type=int, help="Optional limit for number of assemblies to scan")
    parser.add_argument("--min-coverage", type=float, default=80.0)
    parser.add_argument("--min-identity", type=float, default=90.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    legacy_seq = extract_origin_sequence(args.legacy_gb)
    updated_seq = extract_origin_sequence(args.updated_gb)

    legacy_query = write_temp_fasta(legacy_seq, "legacy_operon")
    updated_query = write_temp_fasta(updated_seq, "updated_operon")

    assemblies = sorted(p for p in args.assemblies_dir.glob("*.fasta*"))
    if args.max_assemblies:
        assemblies = assemblies[: args.max_assemblies]

    results: List[MappingResult] = []
    for idx, assembly in enumerate(assemblies, 1):
        print(f"[{idx}/{len(assemblies)}] Processing {assembly.name}")
        legacy_hits = run_blast(legacy_query, assembly, args.threads)
        updated_hits = run_blast(updated_query, assembly, args.threads)
        legacy_hit = next((h for h in legacy_hits if h.coverage >= args.min_coverage and h.pident >= args.min_identity), None)
        updated_hit = next((h for h in updated_hits if h.coverage >= args.min_coverage and h.pident >= args.min_identity), None)
        results.append(MappingResult(assembly=assembly, legacy_hit=legacy_hit, updated_hit=updated_hit))

    # cleanup temp queries
    legacy_query.unlink(missing_ok=True)
    updated_query.unlink(missing_ok=True)

    # Write per-assembly table
    rows = []
    for result in results:
        legacy_stats = alignment_stats(result.legacy_hit)
        updated_stats = alignment_stats(result.updated_hit)
        row = {
            "assembly": result.assembly.name,
            "legacy_found": bool(result.legacy_hit),
            "updated_found": bool(result.updated_hit),
            "legacy_pident": result.legacy_hit.pident if result.legacy_hit else None,
            "legacy_coverage": result.legacy_hit.coverage if result.legacy_hit else None,
            "updated_pident": result.updated_hit.pident if result.updated_hit else None,
            "updated_coverage": result.updated_hit.coverage if result.updated_hit else None,
            "legacy_mismatch_count": legacy_stats["mismatch_count"],
            "legacy_insertion_count": legacy_stats["insertion_count"],
            "legacy_deletion_count": legacy_stats["deletion_count"],
            "legacy_mismatch_positions": legacy_stats["mismatch_positions"],
            "updated_mismatch_count": updated_stats["mismatch_count"],
            "updated_insertion_count": updated_stats["insertion_count"],
            "updated_deletion_count": updated_stats["deletion_count"],
            "updated_mismatch_positions": updated_stats["mismatch_positions"],
            "legacy_qseq": result.legacy_hit.qseq if result.legacy_hit else None,
            "legacy_sseq": result.legacy_hit.sseq if result.legacy_hit else None,
            "updated_qseq": result.updated_hit.qseq if result.updated_hit else None,
            "updated_sseq": result.updated_hit.sseq if result.updated_hit else None,
        }
        rows.append(row)
    df = pd.DataFrame(rows)
    table_path = args.output / "full_operon_mapping.tsv"
    df.to_csv(table_path, sep="\t", index=False)

    summary = summarise_results(results)
    summary_path = args.output / "summary.txt"
    with summary_path.open("w") as handle:
        handle.write("Full operon mapping summary\n")
        handle.write("===========================\n")
        handle.write(f"Legacy GB: {args.legacy_gb}\n")
        handle.write(f"Updated GB: {args.updated_gb}\n")
        handle.write(f"Assemblies scanned: {summary['total']}\n")
        handle.write(f"Legacy hits (>= {args.min_identity}% id, >= {args.min_coverage}% cov): {summary['legacy_hits']}\n")
        handle.write(f"Updated hits (>= {args.min_identity}% id, >= {args.min_coverage}% cov): {summary['updated_hits']}\n")
        handle.write(f"Legacy high-quality hits (>=95/95): {summary['legacy_hq']}\n")
        handle.write(f"Updated high-quality hits (>=95/95): {summary['updated_hq']}\n")

    print(f"Per-assembly table written to {table_path}")
    print(f"Summary written to {summary_path}")


if __name__ == "__main__":
    main()
