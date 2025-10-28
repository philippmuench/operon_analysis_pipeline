#!/usr/bin/env python3
"""Aggregate position-wise mismatch/gap statistics from raw BLAST alignments."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import pandas as pd
from Bio import SeqIO


@dataclass
class BlastAlignment:
    pident: float
    length: int
    bitscore: float
    qstart: int
    qend: int
    qlen: int
    qseq: str
    sseq: str

    @property
    def coverage(self) -> float:
        if self.qlen == 0:
            return 0.0
        return (self.length / self.qlen) * 100.0


def init_counts(length: int) -> Dict[str, List[int]]:
    return {
        "coverage": [0] * (length + 1),
        "mismatch": [0] * (length + 1),
        "deletion": [0] * (length + 1),
        "insertion": [0] * (length + 1),
    }


def ensure_counts(counts: Optional[Dict[str, List[int]]], length: int) -> Dict[str, List[int]]:
    if counts is None:
        return init_counts(length)
    needed = length + 1 - len(counts["coverage"])
    if needed > 0:
        for key in counts:
            counts[key].extend([0] * needed)
    return counts


def accumulate_counts(hit: BlastAlignment, counts: Dict[str, List[int]]):
    qpos = hit.qstart
    step = 1 if hit.qend >= hit.qstart else -1
    max_index = len(counts["coverage"]) - 1

    for q_char, s_char in zip(hit.qseq.upper(), hit.sseq.upper()):
        if q_char == "-" and s_char != "-":
            target = qpos - step
            if target < 1:
                target = 1
            elif target > max_index:
                target = max_index
            counts["insertion"][target] += 1
            continue

        if q_char != "-":
            if s_char == "-":
                counts["coverage"][qpos] += 1
                counts["deletion"][qpos] += 1
            else:
                counts["coverage"][qpos] += 1
                if q_char != s_char:
                    counts["mismatch"][qpos] += 1
            qpos += step


def build_position_table(counts: Optional[Dict[str, List[int]]], reference_seq: str) -> pd.DataFrame:
    if counts is None:
        return pd.DataFrame()

    length = len(reference_seq)
    data = {
        "position": list(range(1, length + 1)),
        "reference_base": list(reference_seq.upper()[:length]),
        "coverage": counts["coverage"][1 : length + 1],
        "mismatch_count": counts["mismatch"][1 : length + 1],
        "deletion_count": counts["deletion"][1 : length + 1],
        "insertion_count": counts["insertion"][1 : length + 1],
    }

    df = pd.DataFrame(data)
    df["mismatch_rate"] = df.apply(
        lambda row: row["mismatch_count"] / row["coverage"] if row["coverage"] > 0 else 0.0,
        axis=1,
    )
    df["deletion_rate"] = df.apply(
        lambda row: row["deletion_count"] / row["coverage"] if row["coverage"] > 0 else 0.0,
        axis=1,
    )
    df["insertion_rate"] = df.apply(
        lambda row: row["insertion_count"] / row["coverage"] if row["coverage"] > 0 else 0.0,
        axis=1,
    )
    return df


def parse_alignment(line: str) -> BlastAlignment:
    parts = line.strip().split("\t")
    return BlastAlignment(
        pident=float(parts[2]),
        length=int(parts[3]),
        bitscore=float(parts[11]),
        qstart=int(parts[6]),
        qend=int(parts[7]),
        qlen=int(parts[12]),
        qseq=parts[13],
        sseq=parts[14],
    )


def load_best_alignment(path: Path, min_identity: float, min_coverage: float) -> Optional[BlastAlignment]:
    best: Optional[BlastAlignment] = None
    for line in path.read_text().splitlines():
        if not line:
            continue
        hit = parse_alignment(line)
        if hit.pident < min_identity:
            continue
        if hit.coverage < min_coverage:
            continue
        if best is None or hit.bitscore > best.bitscore:
            best = hit
    return best


def process_raw_dir(
    raw_dir: Path,
    reference_seq: str,
    min_identity: float,
    min_coverage: float,
) -> pd.DataFrame:
    counts: Optional[Dict[str, List[int]]] = None
    if not raw_dir.exists():
        return pd.DataFrame()

    files = sorted(raw_dir.glob("*.tsv"))
    for path in files:
        hit = load_best_alignment(path, min_identity, min_coverage)
        if hit:
            counts = ensure_counts(counts, hit.qlen)
            accumulate_counts(hit, counts)
    return build_position_table(counts, reference_seq)


def combine_position_tables(legacy: pd.DataFrame, updated: pd.DataFrame) -> pd.DataFrame:
    if legacy.empty and updated.empty:
        return pd.DataFrame()
    merged = legacy.merge(
        updated,
        on="position",
        how="outer",
        suffixes=("_legacy", "_updated"),
    )
    if "reference_base_legacy" in merged.columns:
        merged["reference_base"] = merged["reference_base_legacy"].fillna(merged.get("reference_base_updated"))
        merged.drop(columns=[col for col in merged.columns if col.startswith("reference_base_")], inplace=True)
    else:
        merged.rename(columns={"reference_base_updated": "reference_base"}, inplace=True)
    return merged


def summarise_counts(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()
    return pd.DataFrame(
        {
            "metric": ["coverage", "mismatch_rate", "deletion_rate", "insertion_rate"],
            "mean": [
                df["coverage"].mean(),
                df["mismatch_rate"].mean(),
                df["deletion_rate"].mean(),
                df["insertion_rate"].mean(),
            ],
            "median": [
                df["coverage"].median(),
                df["mismatch_rate"].median(),
                df["deletion_rate"].median(),
                df["insertion_rate"].median(),
            ],
        }
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-raw", type=Path, required=True, help="Path to raw BLAST directory (legacy)")
    parser.add_argument("--updated-raw", type=Path, required=True, help="Path to raw BLAST directory (updated)")
    parser.add_argument("--legacy-reference", type=Path, required=True, help="Legacy reference FASTA/GB containing operon sequence")
    parser.add_argument("--updated-reference", type=Path, required=True, help="Updated reference FASTA/GB containing operon sequence")
    parser.add_argument("--min-identity", type=float, default=90.0)
    parser.add_argument("--min-coverage", type=float, default=80.0)
    parser.add_argument("--output", type=Path, default=Path("position_stats"))
    return parser.parse_args()


def load_reference_sequence(path: Path) -> str:
    if path.suffix.lower() in {".gb", ".gbk"}:
        record = next(SeqIO.parse(str(path), "genbank"))
        return str(record.seq)
    if path.suffix.lower() in {".fa", ".fasta", ".fna"}:
        record = next(SeqIO.parse(str(path), "fasta"))
        return str(record.seq)
    raise ValueError(f"Unsupported reference format: {path}")


def main() -> None:
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    legacy_seq = load_reference_sequence(args.legacy_reference)
    updated_seq = load_reference_sequence(args.updated_reference)

    legacy_df = process_raw_dir(args.legacy_raw, legacy_seq, args.min_identity, args.min_coverage)
    if not legacy_df.empty:
        legacy_df.to_csv(args.output / "legacy_position_summary.tsv", sep="\t", index=False)

    updated_df = process_raw_dir(args.updated_raw, updated_seq, args.min_identity, args.min_coverage)
    if not updated_df.empty:
        updated_df.to_csv(args.output / "updated_position_summary.tsv", sep="\t", index=False)

    combined = combine_position_tables(legacy_df, updated_df)
    if not combined.empty:
        combined.to_csv(args.output / "combined_position_summary.tsv", sep="\t", index=False)

    summary_tables = []
    if not legacy_df.empty:
        summary = summarise_counts(legacy_df)
        summary.insert(0, "dataset", "legacy")
        summary_tables.append(summary)
    if not updated_df.empty:
        summary = summarise_counts(updated_df)
        summary.insert(0, "dataset", "updated")
        summary_tables.append(summary)
    if summary_tables:
        pd.concat(summary_tables, ignore_index=True).to_csv(
            args.output / "position_metric_summary.tsv", sep="\t", index=False
        )


if __name__ == "__main__":
    main()
