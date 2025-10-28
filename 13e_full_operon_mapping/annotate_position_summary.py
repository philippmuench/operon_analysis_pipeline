#!/usr/bin/env python3
"""Annotate operon position summaries with feature metadata and known SNPs."""
from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd


def load_features(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = {"start", "end", "strand", "gene", "locus_tag", "type", "product"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Annotation table {path} missing columns: {sorted(missing)}")
    df = df.copy()
    df.sort_values("start", inplace=True)
    return df



def map_variations(path: Optional[Path]) -> Dict[int, List[str]]:
    if path is None:
        return {}
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path, sep="\t")
    if df.empty or "position" not in df.columns:
        return {}
    mapping: Dict[int, List[str]] = {}
    for _, row in df.iterrows():
        label = str(row.get("label", ""))
        pos_field = str(row["position"])
        numbers = [int(value) for value in re.findall(r"\d+", pos_field)]
        if not numbers:
            continue
        if len(numbers) == 1:
            span = [numbers[0]]
        else:
            span = list(range(numbers[0], numbers[-1] + 1))
        for pos in span:
            mapping.setdefault(pos, []).append(label)
    return mapping


def annotate_positions(
    df: pd.DataFrame,
    features: pd.DataFrame,
    variations: Dict[int, List[str]],
) -> pd.DataFrame:
    annotated = df.copy()

    feature_rows: List[Optional[pd.Series]] = []
    for position in annotated["position"]:
        matched = None
        for _, feature in features.iterrows():
            if int(feature["start"]) <= int(position) <= int(feature["end"]):
                matched = feature
                break
        feature_rows.append(matched)

    def extract(column: str) -> List[Optional[object]]:
        values: List[Optional[object]] = []
        for feature in feature_rows:
            values.append(feature[column] if feature is not None else np.nan)
        return values

    annotated["feature_gene"] = extract("gene")
    annotated["feature_locus_tag"] = extract("locus_tag")
    annotated["feature_product"] = extract("product")
    annotated["feature_type"] = extract("type")
    annotated["feature_start"] = extract("start")
    annotated["feature_end"] = extract("end")
    annotated["feature_strand"] = extract("strand")

    pos_in_feature: List[Optional[int]] = []
    codon_index: List[Optional[int]] = []
    for position, feature_type, strand, start, end in zip(
        annotated["position"],
        annotated["feature_type"],
        annotated["feature_strand"],
        annotated["feature_start"],
        annotated["feature_end"],
    ):
        if pd.isna(start) or pd.isna(end):
            pos_in_feature.append(np.nan)
            codon_index.append(np.nan)
            continue
        if strand == -1:
            offset = int(end) - int(position) + 1
        else:
            offset = int(position) - int(start) + 1
        pos_in_feature.append(offset)
        if feature_type == "CDS" and offset > 0:
            codon = math.ceil(offset / 3.0)
            codon_index.append(codon)
        else:
            codon_index.append(np.nan)

    annotated["feature_position"] = pos_in_feature
    annotated["feature_codon"] = codon_index

    variation_labels: List[str] = []
    for position in annotated["position"]:
        labels = variations.get(int(position))
        variation_labels.append(";".join(labels) if labels else "")
    annotated["known_variation"] = variation_labels

    return annotated


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--combined",
        type=Path,
        default=Path("output_full_run/position_stats/combined_position_summary.tsv"),
        help="Combined position summary to annotate",
    )
    parser.add_argument(
        "--annotation-table",
        type=Path,
        default=Path("../13a_new_reference_operon_extraction/output/operon_genes.tsv"),
        help="Operon feature table with start/end coordinates",
    )
    parser.add_argument(
        "--variation-table",
        type=Path,
        default=Path("../13a_new_reference_operon_extraction/output/operon_variations.tsv"),
        help="Optional table listing known variation positions",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Path to write annotated table (default: overwrite combined file)",
    )
    args = parser.parse_args()

    combined_path = args.combined
    if not combined_path.exists():
        raise FileNotFoundError(combined_path)

    df = pd.read_csv(combined_path, sep="\t")
    if "position" not in df.columns:
        raise ValueError(f"Expected 'position' column in {combined_path}")

    features = load_features(args.annotation_table)
    variations = map_variations(args.variation_table if args.variation_table.exists() else None)

    annotated = annotate_positions(df, features, variations)

    output_path = args.output or combined_path
    annotated.to_csv(output_path, sep="\t", index=False)


if __name__ == "__main__":
    main()
