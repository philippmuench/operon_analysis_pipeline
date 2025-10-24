#!/usr/bin/env python3
"""Summarise allele usage at predefined operon positions across BLAST outputs."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

BASES = {"A", "C", "G", "T", "N"}
INSERT_ABSENT = "absent"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate allele calls for selected operon positions from BLAST output (.tsv)."
    )
    parser.add_argument(
        "raw_dir",
        type=Path,
        help="Directory containing BLAST tabular outputs (default: updated mapping dir).",
        nargs="?",
        default=Path("13e_full_operon_mapping/output_full_run/raw_blast/updated"),
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=0.0,
        help="Minimum percent identity for an alignment to be considered (default: 0).",
    )
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=0.0,
        help="Minimum query coverage (percent) per subject to be retained (default: 0).",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional limit on number of TSV files to process (for testing).",
    )
    parser.add_argument(
        "--site-table",
        type=Path,
        help="Optional path to write per-site allele counts (CSV if *.csv, else TSV).",
    )
    parser.add_argument(
        "--site-calls",
        type=Path,
        help="Optional path to write per-assembly allele calls per site (CSV if *.csv, else TSV).",
    )
    parser.add_argument(
        "--site-summary",
        type=Path,
        help="Optional path to write per-site summary statistics (CSV if *.csv, else TSV).",
    )
    parser.add_argument(
        "--figure-out",
        type=Path,
        help="Optional path to save a stacked bar chart of allele usage (requires matplotlib).",
    )
    return parser.parse_args()


@dataclass(frozen=True)
class VariantSite:
    label: str
    kind: str  # snv, block, insertion
    positions: Tuple[int, ...]
    origin: str
    variant: str


VARIANT_SITES: Tuple[VariantSite, ...] = (
    VariantSite("859", "snv", (859,), "G", "A"),
    VariantSite("2155", "snv", (2155,), "C", "A"),
    VariantSite("3455", "snv", (3455,), "G", "T"),
    VariantSite("3568", "snv", (3568,), "C", "G"),
    VariantSite("4652..4653", "block", (4652, 4653), "CT", "TG"),
    VariantSite("4659", "snv", (4659,), "T", "A"),
    VariantSite("4673", "snv", (4673,), "C", "G"),
    VariantSite("7501^7502", "insertion", (7501, 7502), INSERT_ABSENT, "A"),
)

TARGET_BASE_POSITIONS = {pos for site in VARIANT_SITES for pos in site.positions if site.kind != "insertion"}
TARGET_INSERT_KEYS = {tuple(sorted(site.positions)) for site in VARIANT_SITES if site.kind == "insertion"}
SITE_BY_LABEL = {site.label: site for site in VARIANT_SITES}


@dataclass
class Alignment:
    subject: str
    pident: float
    length: int
    mismatch: int
    gaps: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    qlen: int
    qseq: str
    sseq: str

    @property
    def covered_bases(self) -> int:
        return sum(1 for base in self.qseq if base != "-")


@dataclass
class SubjectBundle:
    alignments: List[Alignment]
    qlen: int
    total_covered: int
    total_bitscore: float


def parse_alignment(line: str) -> Alignment:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 15:
        raise ValueError("Unexpected alignment format; expected at least 15 columns.")
    return Alignment(
        subject=fields[1],
        pident=float(fields[2]),
        length=int(fields[3]),
        mismatch=int(fields[4]),
        gaps=int(fields[5]),
        qstart=int(fields[6]),
        qend=int(fields[7]),
        sstart=int(fields[8]),
        send=int(fields[9]),
        evalue=float(fields[10]),
        bitscore=float(fields[11]),
        qlen=int(fields[12]),
        qseq=fields[13].upper(),
        sseq=fields[14].upper(),
    )


def load_subject_bundles(path: Path, min_identity: float) -> Dict[str, SubjectBundle]:
    bundles: Dict[str, SubjectBundle] = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            aln = parse_alignment(line)
            if aln.pident < min_identity:
                continue
            bundle = bundles.get(aln.subject)
            if bundle is None:
                bundle = SubjectBundle([], aln.qlen, 0, 0.0)
                bundles[aln.subject] = bundle
            bundle.alignments.append(aln)
            bundle.total_covered += aln.covered_bases
            bundle.total_bitscore += aln.bitscore
    return bundles


def select_best_bundle(bundles: Dict[str, SubjectBundle], min_coverage: float) -> SubjectBundle | None:
    best: SubjectBundle | None = None
    for bundle in bundles.values():
        if bundle.qlen <= 0:
            continue
        coverage_pct = 100.0 * (bundle.total_covered / bundle.qlen)
        if coverage_pct < min_coverage:
            continue
        if best is None:
            best = bundle
            continue
        if bundle.total_covered > best.total_covered:
            best = bundle
            continue
        if bundle.total_covered == best.total_covered and bundle.total_bitscore > best.total_bitscore:
            best = bundle
    return best


def summarise_base(values: Sequence[str]) -> str:
    if not values:
        return "missing"
    alleles = {val for val in values if val}
    if not alleles:
        return "missing"
    if alleles == {"-"}:
        return "gap"
    if len(alleles) == 1:
        return alleles.pop()
    return "ambiguous:" + "/".join(sorted(alleles))


def summarise_block(bases: Sequence[str]) -> str:
    if "missing" in bases:
        return "missing"
    if "gap" in bases:
        return "gap"
    clean: List[str] = []
    for base in bases:
        if base.startswith("ambiguous"):
            return "ambiguous"
        if base not in BASES:
            return "ambiguous"
        clean.append(base)
    return "".join(clean)


def summarise_insertion(seqs: Iterable[str]) -> str:
    collected = [seq for seq in seqs if seq]
    if not collected:
        return INSERT_ABSENT
    unique = sorted({seq.upper() for seq in collected})
    if len(unique) == 1:
        return unique[0]
    return "ambiguous:" + "/".join(unique)


def extract_calls(bundle: SubjectBundle) -> Tuple[Dict[int, str], Dict[Tuple[int, int], List[str]]]:
    base_values: Dict[int, List[str]] = {pos: [] for pos in TARGET_BASE_POSITIONS}
    insertion_values: Dict[Tuple[int, int], List[str]] = defaultdict(list)

    for aln in bundle.alignments:
        q_pos = aln.qstart
        q_step = 1 if aln.qend >= aln.qstart else -1
        prev_q_pos: int | None = None
        current_ins_key: Tuple[int, int] | None = None
        current_ins_seq: List[str] = []

        def flush_insertion():
            nonlocal current_ins_key, current_ins_seq
            if current_ins_key is not None and current_ins_seq:
                insertion_values[current_ins_key].append("".join(current_ins_seq))
            current_ins_key = None
            current_ins_seq = []

        for q_char, s_char in zip(aln.qseq, aln.sseq):
            q_curr: int | None = None
            if q_char != "-":
                q_curr = q_pos
                q_pos += q_step
            if q_char == "-" and s_char != "-" and prev_q_pos is not None:
                key = tuple(sorted((prev_q_pos, q_pos)))
                if key in TARGET_INSERT_KEYS:
                    if current_ins_key != key:
                        flush_insertion()
                        current_ins_key = key
                    current_ins_seq.append(s_char)
                continue
            flush_insertion()
            if q_curr is not None and q_curr in base_values:
                base_values[q_curr].append(s_char if s_char != "-" else "-")
                prev_q_pos = q_curr
            elif q_curr is not None:
                prev_q_pos = q_curr
        flush_insertion()

    base_calls = {pos: summarise_base(values) for pos, values in base_values.items()}
    return base_calls, insertion_values


def aggregate_sites(base_calls: Dict[int, str], insertion_values: Dict[Tuple[int, int], List[str]]) -> Dict[str, str]:
    site_calls: Dict[str, str] = {}
    for site in VARIANT_SITES:
        if site.kind == "snv":
            call = base_calls.get(site.positions[0], "missing")
        elif site.kind == "block":
            call = summarise_block([base_calls.get(pos, "missing") for pos in site.positions])
        elif site.kind == "insertion":
            key = tuple(sorted(site.positions))
            call = summarise_insertion(insertion_values.get(key, []))
        else:
            raise ValueError(f"Unsupported site kind: {site.kind}")
        site_calls[site.label] = call
    return site_calls


def prepare_summary_tables(counters: Dict[str, Counter]) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    site_summaries: List[Dict[str, object]] = []
    site_details: List[Dict[str, object]] = []

    for site in VARIANT_SITES:
        counter = counters[site.label]
        total = sum(counter.values())
        origin_count = counter.get(site.origin, 0)
        variant_count = counter.get(site.variant, 0)
        gap_count = counter.get("gap", 0)
        missing_count = counter.get("missing", 0)
        other_items = [(allele, count) for allele, count in counter.items() if allele not in {site.origin, site.variant, "gap", "missing"}]
        other_items.sort()
        other_total = sum(count for _, count in other_items)

        def pct(count: int) -> float:
            return round((count / total * 100.0), 4) if total else 0.0

        other_breakdown = ";".join(
            f"{allele}:{count} ({pct(count):.2f}%)" for allele, count in other_items
        )

        site_summaries.append(
            {
                "site": site.label,
                "origin": site.origin,
                "variant": site.variant,
                "total_calls": total,
                "origin_count": origin_count,
                "origin_percent": pct(origin_count),
                "variant_count": variant_count,
                "variant_percent": pct(variant_count),
                "gap_count": gap_count,
                "gap_percent": pct(gap_count),
                "missing_count": missing_count,
                "missing_percent": pct(missing_count),
                "other_count": other_total,
                "other_percent": pct(other_total),
                "other_breakdown": other_breakdown,
            }
        )

        recorded: Dict[str, bool] = {}
        for allele, count in sorted(counter.items()):
            if allele == site.origin:
                category = "origin"
            elif allele == site.variant:
                category = "variant"
            elif allele == "missing":
                category = "missing"
            elif allele == "gap":
                category = "gap"
            else:
                category = "other"
            site_details.append(
                {
                    "site": site.label,
                    "allele": allele,
                    "count": count,
                    "percent": pct(count),
                    "category": category,
                    "is_origin": allele == site.origin,
                    "is_variant": allele == site.variant,
                }
            )
            recorded[allele] = True
        for allele, is_origin in ((site.origin, True), (site.variant, False)):
            if allele and allele not in recorded:
                site_details.append(
                    {
                        "site": site.label,
                        "allele": allele,
                        "count": 0,
                        "percent": 0.0,
                        "category": "origin" if is_origin else "variant",
                        "is_origin": is_origin,
                        "is_variant": (not is_origin),
                    }
                )

    return site_summaries, site_details


def process_directory(
    raw_dir: Path,
    min_identity: float,
    min_coverage: float,
    limit: int | None,
) -> Tuple[Dict[str, Counter], int, List[Dict[str, object]]]:
    counters: Dict[str, Counter] = {site.label: Counter() for site in VARIANT_SITES}
    files = sorted(raw_dir.glob("*.tsv"))
    processed = 0
    call_records: List[Dict[str, object]] = []

    for idx, path in enumerate(files, start=1):
        if limit is not None and idx > limit:
            break
        bundles = load_subject_bundles(path, min_identity)
        if not bundles:
            continue
        bundle = select_best_bundle(bundles, min_coverage)
        if bundle is None:
            continue
        base_calls, insertion_values = extract_calls(bundle)
        site_calls = aggregate_sites(base_calls, insertion_values)
        for label, call in site_calls.items():
            counters[label][call] += 1
            site = SITE_BY_LABEL[label]
            if call == site.origin:
                category = "origin"
            elif call == site.variant:
                category = "variant"
            elif call == "missing":
                category = "missing"
            elif call == "gap":
                category = "gap"
            else:
                category = "other"
            call_records.append(
                {
                    "assembly": path.stem,
                    "site": label,
                    "allele": call,
                    "category": category,
                    "is_origin": call == site.origin,
                    "is_variant": call == site.variant,
                }
            )
        processed += 1
    return counters, processed, call_records


def render_summary(counters: Dict[str, Counter], processed: int) -> str:
    lines = [f"Processed assemblies: {processed}"]
    for site in VARIANT_SITES:
        counter = counters[site.label]
        total = sum(counter.values())
        origin_count = counter.get(site.origin, 0)
        variant_count = counter.get(site.variant, 0)
        gap_count = counter.get("gap", 0)
        missing_count = counter.get("missing", 0)
        others = {
            allele: count
            for allele, count in counter.items()
            if allele not in {site.origin, site.variant, "gap", "missing"}
        }
        lines.append(
            f"\nSite {site.label} (origin={site.origin}, variant={site.variant})\n"
            f"  calls counted : {total}\n"
            f"  match origin  : {origin_count}\n"
            f"  match variant : {variant_count}\n"
            f"  gaps          : {gap_count}\n"
            f"  missing       : {missing_count}"
        )
        if others:
            formatted = ", ".join(f"{allele}={count}" for allele, count in sorted(others.items()))
            lines.append(f"  others        : {formatted}")
    return "\n".join(lines)


def write_table(path: Path, fieldnames: Sequence[str], records: Sequence[Dict[str, object]]) -> None:
    delimiter = "," if path.suffix.lower() == ".csv" else "\t"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter=delimiter,
            lineterminator="\n",
        )
        writer.writeheader()
        for record in records:
            writer.writerow(record)


def plot_counts(path: Path, summaries: Sequence[Dict[str, object]]) -> None:
    if not summaries:
        # Still create an empty file to signal execution, but no plotting needed.
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
        return

    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - runtime dependency
        raise SystemExit("matplotlib is required to produce figures; install it or omit --figure-out.") from exc

    categories = ["origin", "variant", "other", "gap", "missing"]
    colors = {
        "origin": "#1f77b4",
        "variant": "#d62728",
        "other": "#ff7f0e",
        "gap": "#9467bd",
        "missing": "#7f7f7f",
    }

    labels: List[str] = []
    values: Dict[str, List[int]] = {key: [] for key in categories}

    for summary in summaries:
        labels.append(str(summary["site"]))
        values["origin"].append(int(summary["origin_count"]))
        values["variant"].append(int(summary["variant_count"]))
        values["other"].append(int(summary["other_count"]))
        values["gap"].append(int(summary["gap_count"]))
        values["missing"].append(int(summary["missing_count"]))

    path.parent.mkdir(parents=True, exist_ok=True)
    fig_width = max(6.0, 1.2 * len(labels))
    fig, ax = plt.subplots(figsize=(fig_width, 4.5))

    bottom = [0] * len(labels)
    for category in categories:
        heights = values[category]
        if not any(heights):
            continue
        ax.bar(labels, heights, bottom=bottom, label=category, color=colors.get(category, None))
        bottom = [b + h for b, h in zip(bottom, heights)]

    ax.set_xlabel("Variant site")
    ax.set_ylabel("Assemblies")
    ax.set_title("Allele usage across mapped assemblies")
    ax.legend()
    ax.set_ylim(0, max(bottom) * 1.05 if bottom else 1)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    if not args.raw_dir.exists():
        raise SystemExit(f"Input directory not found: {args.raw_dir}")
    counters, processed, call_records = process_directory(
        raw_dir=args.raw_dir,
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
        limit=args.limit,
    )
    summary_text = render_summary(counters, processed)
    print(summary_text)

    site_summaries, site_details = prepare_summary_tables(counters)

    if args.site_summary:
        write_table(
            args.site_summary,
            (
                "site",
                "origin",
                "variant",
                "total_calls",
                "origin_count",
                "origin_percent",
                "variant_count",
                "variant_percent",
                "gap_count",
                "gap_percent",
                "missing_count",
                "missing_percent",
                "other_count",
                "other_percent",
                "other_breakdown",
            ),
            site_summaries,
        )

    if args.site_table:
        write_table(
            args.site_table,
            ("site", "allele", "count", "percent", "category", "is_origin", "is_variant"),
            site_details,
        )

    if args.site_calls:
        write_table(
            args.site_calls,
            ("assembly", "site", "allele", "category", "is_origin", "is_variant"),
            call_records,
        )

    if args.figure_out:
        plot_counts(args.figure_out, site_summaries)


if __name__ == "__main__":
    main()
