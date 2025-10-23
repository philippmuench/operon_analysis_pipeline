#!/usr/bin/env python3
"""Extract the fructoselysine/glucoselysine operon genes from the SNP-annotated GenBank file.

The script mirrors the behaviour of `02_reference_operon_extraction/extract_operon_genes.py`
but targets the GenBank file located in this directory and keeps all outputs local to
`13_operon/`.  Gene coordinates are based on the curated annotations embedded in
`FL_operon_with_SNPs.gb`.
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path
from typing import Dict, Iterable, List

import pandas as pd
from Bio import SeqIO
from Bio import BiopythonParserWarning
from Bio.Seq import Seq

PREFERRED_NOTE_ORDER = (
    "geneious type: kofam",
    "geneious type: eggnogdb",
    "geneious type: ncbi",
    "geneious type: operonmapper",
    "geneious type: rgi",
)

ANNOTATION_RULES: List[Dict[str, object]] = [
    {
        "keywords": {"fructoselysine-6-phosphate deglycase"},
        "gene": "frpC",
        "product": "fructoselysine-6-phosphate deglycase",
        "feature_type": "CDS",
    },
    {
        "keywords": {"glucoselysine-6-phosphate deglycase"},
        "gene": "glpC",
        "product": "glucoselysine-6-phosphate deglycase",
        "feature_type": "CDS",
    },
    {
        "keywords": {
            "pts system mannose/fructose/sorbose family iid component",
            "fructoselysine/glucoselysine pts system eiid component",
            "pts system mannose/fructose/sorbose family transporter subunit iid",
        },
        "gene": "ptsD",
        "product": "PTS system mannose/fructose/sorbose family IID component",
        "feature_type": "CDS",
    },
    {
        "keywords": {
            "pts system sorbose-specific iic component",
            "fructoselysine/glucoselysine pts system eiic component",
            "pts sugar transporter subunit iic",
        },
        "gene": "ptsC",
        "product": "PTS system sorbose-specific IIC component",
        "feature_type": "CDS",
    },
    {
        "keywords": {
            "pts system sorbose subfamily iib component",
            "fructoselysine/glucoselysine pts system eiib component",
            "pts sugar transporter subunit iib",
        },
        "gene": "ptsB",
        "product": "PTS system sorbose subfamily IIB component",
        "feature_type": "CDS",
    },
    {
        "keywords": {
            "pts system fructose iia component",
            "fructoselysine/glucoselysine pts system eiia component",
            "pts mannose transporter subunit iia",
        },
        "gene": "ptsA",
        "product": "PTS system fructose IIA component",
        "feature_type": "CDS",
    },
    {
        "keywords": {
            "sigma-54 dependent transcriptional regulator",
            "sigma 54-interacting transcriptional regulator",
            "sigma-54 interaction domain",
        },
        "gene": "fruR",
        "product": "sigma-54 dependent transcriptional regulator",
        "feature_type": "CDS",
    },
    {
        "keywords": {"-12 region"},
        "gene": "promoter",
        "product": "-12 region",
        "feature_type": "regulatory",
    },
    {
        "keywords": {"pribnow box"},
        "gene": "pribnow_box",
        "product": "Pribnow box (-10 box)",
        "feature_type": "regulatory",
    },
]


def translate_sequence(sequence: Seq) -> str:
    """Translate a nucleotide sequence while tolerating trailing incomplete codons."""
    try:
        return str(sequence.translate(to_stop=True))
    except Exception:
        return str(sequence.translate())


def normalize_identifier(text: str) -> str:
    """Generate a filesystem-friendly identifier from free text."""
    cleaned = text.lower().replace("/", "_").replace("-", "_").replace(" ", "_")
    cleaned = "".join(ch for ch in cleaned if ch.isalnum() or ch == "_")
    cleaned = cleaned.strip("_")
    return cleaned or "feature"


def choose_annotation_label(annotations: List[Dict[str, List[str]]]) -> str:
    """Pick the most informative label available for a feature."""
    for priority in PREFERRED_NOTE_ORDER:
        for quals in annotations:
            for note in quals.get("note", []):
                if priority in note.lower():
                    names = quals.get("standard_name")
                    if names:
                        return names[0]

    for quals in annotations:
        names = quals.get("standard_name")
        if names:
            return names[0]

    for quals in annotations:
        funcs = quals.get("Function")
        if funcs:
            return funcs[0]

    return ""


def collect_candidate_strings(annotations: List[Dict[str, List[str]]], label: str) -> List[str]:
    candidates: List[str] = []
    if label:
        candidates.append(label)

    for quals in annotations:
        for key in ("standard_name", "Function", "product", "gene", "note"):
            values = quals.get(key)
            if not values:
                continue
            for value in values:
                if value not in candidates:
                    candidates.append(value)

    return candidates


def classify_feature(candidates: Iterable[str], length_nt: int) -> Dict[str, str]:
    """Map raw annotation strings to gene name, product, and feature type."""
    lowered = [candidate.lower() for candidate in candidates if candidate]

    for rule in ANNOTATION_RULES:
        if any(keyword in candidate for candidate in lowered for keyword in rule["keywords"]):
            return {
                "gene": rule["gene"],
                "product": rule["product"],
                "type": rule["feature_type"],
            }

    fallback = next((candidate for candidate in candidates if candidate), "unknown feature")
    feature_type = "regulatory" if length_nt <= 15 else "misc_feature"
    if "box" in fallback.lower() or "region" in fallback.lower():
        feature_type = "regulatory"

    return {
        "gene": normalize_identifier(fallback),
        "product": fallback,
        "type": feature_type,
    }


def collect_misc_features(record) -> List[Dict[str, object]]:
    """Group duplicate misc_feature annotations that share coordinates."""
    grouped: Dict[tuple, Dict[str, object]] = {}

    for feature in record.features:
        if feature.type != "misc_feature":
            continue

        start = int(feature.location.start) + 1
        end = int(feature.location.end)
        key = (start, end)

        bucket = grouped.setdefault(
            key,
            {
                "start": start,
                "end": end,
                "strand": feature.location.strand or 1,
                "annotations": [],
            },
        )
        bucket["annotations"].append(feature.qualifiers)

    return [grouped[key] for key in sorted(grouped)]


def extract_operon_entries(record) -> pd.DataFrame:
    """Build a tidy table with sequence content for each operon element."""
    entries = []
    cds_counter = 0

    for feature in collect_misc_features(record):
        start = feature["start"]
        end = feature["end"]
        strand = feature["strand"]
        length_nt = end - start + 1

        label = choose_annotation_label(feature["annotations"])
        candidates = collect_candidate_strings(feature["annotations"], label)
        classification = classify_feature(candidates, length_nt)

        subseq = record.seq[start - 1 : end]
        if strand == -1:
            subseq = subseq.reverse_complement()

        translation = ""
        locus_tag = f"fl_operon_{classification['gene']}"
        if classification["type"] == "CDS":
            translation = translate_sequence(subseq)
            cds_counter += 1
            locus_tag = f"fl_operon_CDS_{cds_counter}"

        entries.append(
            {
                "locus_tag": locus_tag,
                "gene": classification["gene"],
                "product": classification["product"],
                "start": start,
                "end": end,
                "strand": strand,
                "type": classification["type"],
                "length_nt": length_nt,
                "sequence_nt": str(subseq),
                "translation": translation,
            }
        )

    return pd.DataFrame(entries)


def extract_variations(record) -> pd.DataFrame:
    """Collect SNP annotations (type `variation`) from the GenBank record."""
    variations: List[Dict[str, object]] = []

    for feature in record.features:
        if feature.type != "variation":
            continue

        start = int(feature.location.start) + 1  # Convert to 1-based coordinate.
        end = int(feature.location.end)
        # Single base variations have start == end; represent as string for readability.
        label = feature.qualifiers.get("standard_name", feature.qualifiers.get("note", ["unknown"]))[0]
        variations.append(
            {
                "position": start if start == end else f"{start}-{end}",
                "label": label,
            }
        )

    return pd.DataFrame(variations)


def write_fasta(records: pd.DataFrame, fasta_path: Path, sequence_column: str) -> None:
    """Write FASTA with operon entries from the provided DataFrame."""
    with fasta_path.open("w") as handle:
        for _, row in records.iterrows():
            sequence = row[sequence_column]
            if not sequence:
                continue

            header = f"{row['locus_tag']}|{row['gene']}|{row['product']}".replace(" ", "_").replace(",", "")
            handle.write(f">{header}\n{sequence}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--gb",
        type=Path,
        default=Path(__file__).with_name("FL_operon_with_SNPs.gb"),
        help="GenBank file containing the operon (default: %(default)s)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("output"),
        help="Directory for output tables and FASTA files (default: %(default)s)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    gb_path: Path = args.gb
    output_dir: Path = args.output

    warnings.simplefilter("ignore", BiopythonParserWarning)

    output_dir.mkdir(parents=True, exist_ok=True)

    record = next(SeqIO.parse(str(gb_path), "genbank"))
    entries = extract_operon_entries(record)

    table_path = output_dir / "operon_genes.tsv"
    entries.drop(columns=["sequence_nt"], inplace=False).to_csv(table_path, sep="\t", index=False)

    coding_entries = entries[entries["type"] == "CDS"].copy()
    regulatory_entries = entries[entries["type"] != "CDS"].copy()

    write_fasta(coding_entries, output_dir / "operon_genes_nt.fasta", "sequence_nt")
    write_fasta(coding_entries, output_dir / "operon_genes_protein.fasta", "translation")
    write_fasta(regulatory_entries, output_dir / "operon_regulatory_nt.fasta", "sequence_nt")

    variations = extract_variations(record)
    if not variations.empty:
        variations.to_csv(output_dir / "operon_variations.tsv", sep="\t", index=False)

    print(f"Wrote gene table to {table_path}")
    print("Coding entries:")
    for _, row in coding_entries.iterrows():
        print(f"  {row['gene']}: {row['product']} ({row['length_nt']} nt)")

    if not regulatory_entries.empty:
        print("Regulatory entries:")
        for _, row in regulatory_entries.iterrows():
            print(f"  {row['gene']}: {row['product']} ({row['length_nt']} nt)")

    if not variations.empty:
        print(f"Recorded {len(variations)} variation annotations.")


if __name__ == "__main__":
    main()
