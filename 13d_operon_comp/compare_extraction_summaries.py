#!/usr/bin/env python3
"""Compare operon extraction pipeline summaries from two runs."""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd


@dataclass
class SummaryData:
    path: Path
    generated: Optional[str]
    sequences: Dict[str, int]
    alignments: Dict[str, Tuple[int, Optional[int]]]
    top_scores: List[Tuple[str, float]]


SEQUENCE_PATTERN = re.compile(r"^(?P<gene>[A-Za-z0-9_]+):\s+([0-9,]+)\s+sequences extracted$")
ALIGNMENT_PATTERN = re.compile(
    r"^(?P<gene>[A-Za-z0-9_]+):\s+([0-9,]+)\s+sequences,\s+([0-9,]+)$"
)
ALIGNMENT_ALT_PATTERN = re.compile(r"^(?P<gene>[A-Za-z0-9_]+):\s+([0-9,]+)\s+sequences aligned$")
ALIGNMENT_LENGTH_PATTERN = re.compile(r"^([0-9,]+)\s+bp alignment$")
TOP_SCORE_PATTERN = re.compile(r"^(?P<gene>[A-Za-z0-9_]+):\s+([0-9.]+)")


def parse_summary(path: Path) -> SummaryData:
    lines = path.read_text().splitlines()
    generated: Optional[str] = None
    sequences: Dict[str, int] = {}
    alignments: Dict[str, Tuple[int, Optional[int]]] = {}
    top_scores: List[Tuple[str, float]] = []

    section: Optional[str] = None
    for idx, raw in enumerate(lines):
        line = raw.strip()
        if not line:
            continue

        if line.startswith("Generated:"):
            generated = line.split("Generated:", 1)[1].strip()
            continue

        if line.startswith("PRIMARY ANALYSIS"):
            section = "primary"
            continue
        if line.startswith("ALIGNMENTS"):
            section = "alignments"
            continue
        if line.startswith("CONSERVATION METRICS"):
            section = "metrics"
            continue

        if section == "primary":
            match = SEQUENCE_PATTERN.match(line)
            if match:
                gene = match.group("gene").lower()
                count = int(match.group(0).split(":")[1].split()[0].replace(",", ""))
                sequences[gene] = count
        elif section == "alignments":
            match = ALIGNMENT_PATTERN.match(line)
            if match:
                gene = match.group("gene").lower()
                seq_count = int(match.group(2).replace(",", ""))
                # Alignment length may be in same line but usually next line
                alignment_length = int(match.group(3).replace(",", ""))
                next_line = lines[idx + 1].strip() if idx + 1 < len(lines) else ""
                match_length = ALIGNMENT_LENGTH_PATTERN.match(next_line)
                if match_length:
                    alignment_length = int(match_length.group(1).replace(",", ""))
                alignments[gene] = (seq_count, alignment_length)
            else:
                match_alt = ALIGNMENT_ALT_PATTERN.match(line)
                if match_alt:
                    gene = match_alt.group("gene").lower()
                    seq_count = int(match_alt.group(2).replace(",", ""))
                    alignments[gene] = (seq_count, None)
        elif section == "metrics":
            if line.startswith("Top conservation scores"):
                continue
            match = TOP_SCORE_PATTERN.match(line)
            if match:
                gene = match.group("gene").lower()
                score = float(match.group(2))
                top_scores.append((gene, score))

    return SummaryData(path=path, generated=generated, sequences=sequences, alignments=alignments, top_scores=top_scores)


def build_sequences_table(a: SummaryData, b: SummaryData) -> pd.DataFrame:
    genes = sorted(set(a.sequences.keys()) | set(b.sequences.keys()))
    rows = []
    for gene in genes:
        count_a = a.sequences.get(gene)
        count_b = b.sequences.get(gene)
        diff = None
        if count_a is not None and count_b is not None:
            diff = count_b - count_a
        rows.append({"gene": gene, "run_a_sequences": count_a, "run_b_sequences": count_b, "difference": diff})
    return pd.DataFrame(rows)


def build_alignments_table(a: SummaryData, b: SummaryData) -> pd.DataFrame:
    genes = sorted(set(a.alignments.keys()) | set(b.alignments.keys()))
    rows = []
    for gene in genes:
        seq_a, len_a = a.alignments.get(gene, (None, None))
        seq_b, len_b = b.alignments.get(gene, (None, None))
        diff_seq = (seq_b - seq_a) if None not in (seq_a, seq_b) else None
        diff_len = (len_b - len_a) if None not in (len_a, len_b) else None
        rows.append(
            {
                "gene": gene,
                "run_a_sequences": seq_a,
                "run_b_sequences": seq_b,
                "sequence_difference": diff_seq,
                "run_a_alignment_len": len_a,
                "run_b_alignment_len": len_b,
                "alignment_len_difference": diff_len,
            }
        )
    return pd.DataFrame(rows)


def build_top_scores_table(a: SummaryData, b: SummaryData) -> pd.DataFrame:
    def to_df(summary: SummaryData, label: str) -> pd.DataFrame:
        return pd.DataFrame(summary.top_scores, columns=["gene", label]).set_index("gene")

    df_a = to_df(a, "run_a_score") if a.top_scores else pd.DataFrame()
    df_b = to_df(b, "run_b_score") if b.top_scores else pd.DataFrame()
    if df_a.empty and df_b.empty:
        return pd.DataFrame()
    merged = df_a.join(df_b, how="outer")
    merged["difference"] = merged.apply(
        lambda row: row["run_b_score"] - row["run_a_score"] if pd.notna(row["run_a_score"]) and pd.notna(row["run_b_score"]) else pd.NA,
        axis=1,
    )
    return merged.reset_index()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary_a", type=Path, help="Path to first extraction summary")
    parser.add_argument("summary_b", type=Path, help="Path to second extraction summary")
    parser.add_argument("--labels", nargs=2, metavar=("LABEL_A", "LABEL_B"), help="Labels for the two runs")
    parser.add_argument("--output-dir", type=Path, default=Path("output"), help="Directory for comparison outputs")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    summary_a = parse_summary(args.summary_a)
    summary_b = parse_summary(args.summary_b)

    label_a, label_b = args.labels if args.labels else ("run_a", "run_b")

    seq_df = build_sequences_table(summary_a, summary_b)
    seq_df.rename(columns={"run_a_sequences": f"{label_a}_sequences", "run_b_sequences": f"{label_b}_sequences"}, inplace=True)

    align_df = build_alignments_table(summary_a, summary_b)
    align_df.rename(
        columns={
            "run_a_sequences": f"{label_a}_sequences",
            "run_b_sequences": f"{label_b}_sequences",
            "run_a_alignment_len": f"{label_a}_alignment_len",
            "run_b_alignment_len": f"{label_b}_alignment_len",
        },
        inplace=True,
    )

    top_df = build_top_scores_table(summary_a, summary_b)
    if not top_df.empty:
        top_df.rename(columns={"run_a_score": f"{label_a}_score", "run_b_score": f"{label_b}_score"}, inplace=True)

    seq_path = args.output_dir / "sequences_comparison.tsv"
    seq_df.to_csv(seq_path, sep="\t", index=False)

    align_path = args.output_dir / "alignments_comparison.tsv"
    align_df.to_csv(align_path, sep="\t", index=False)

    if not top_df.empty:
        top_path = args.output_dir / "top_conservation_comparison.tsv"
        top_df.to_csv(top_path, sep="\t", index=False)

    summary_path = args.output_dir / "summary.txt"
    with summary_path.open("w") as handle:
        handle.write("Operon extraction comparison\n")
        handle.write("============================\n\n")
        handle.write(f"Run A ({label_a}): {summary_a.path}\n")
        if summary_a.generated:
            handle.write(f"  Generated: {summary_a.generated}\n")
        handle.write(f"Run B ({label_b}): {summary_b.path}\n")
        if summary_b.generated:
            handle.write(f"  Generated: {summary_b.generated}\n")
        handle.write("\nSequences per gene written to: \n  " + str(seq_path) + "\n")
        handle.write("Alignments comparison written to:\n  " + str(align_path) + "\n")
        if not top_df.empty:
            handle.write("Top conservation scores table: \n  " + str(top_path) + "\n")

    print(f"Saved sequences comparison to {seq_path}")
    print(f"Saved alignments comparison to {align_path}")
    if not top_df.empty:
        print(f"Saved top conservation comparison to {top_path}")
    print(f"Summary written to {summary_path}")


if __name__ == "__main__":
    main()
