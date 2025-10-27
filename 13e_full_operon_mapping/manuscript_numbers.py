#!/usr/bin/env python3
"""Generate reproducible summary statistics for the full operon mapping run."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from statistics import mean
from typing import Iterable, List


@dataclass
class MappingStats:
    total: int
    legacy_hits: int
    updated_hits: int
    both_hits: int
    legacy_only: int
    updated_only: int
    no_hits: int
    legacy_hq: int
    updated_hq: int
    identity_deltas: List[float]
    coverage_deltas: List[float]


@dataclass
class PositionMetric:
    position: int
    reference_base: str
    coverage: int
    mismatch_rate: float


@dataclass
class VariantSite:
    label: str
    origin: str
    variant: str
    total_calls: int
    variant_count: int
    variant_percent: float


def parse_bool(value: str) -> bool:
    return value.strip().lower() in {"true", "1", "yes"}


def parse_float(value: str) -> float | None:
    value = value.strip()
    if not value or value.lower() in {"na", "nan"}:
        return None
    return float(value)


def load_mapping_stats(path: Path) -> MappingStats:
    total = legacy_hits = updated_hits = both_hits = 0
    legacy_only = updated_only = no_hits = 0
    legacy_hq = updated_hq = 0
    identity_deltas: List[float] = []
    coverage_deltas: List[float] = []

    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            total += 1
            legacy_found = parse_bool(row["legacy_found"])
            updated_found = parse_bool(row["updated_found"])
            legacy_hits += int(legacy_found)
            updated_hits += int(updated_found)

            if legacy_found and updated_found:
                both_hits += 1
            elif legacy_found:
                legacy_only += 1
            elif updated_found:
                updated_only += 1
            else:
                no_hits += 1

            legacy_pident = parse_float(row.get("legacy_pident", ""))
            legacy_cov = parse_float(row.get("legacy_coverage", ""))
            updated_pident = parse_float(row.get("updated_pident", ""))
            updated_cov = parse_float(row.get("updated_coverage", ""))

            if legacy_found and legacy_pident is not None and legacy_cov is not None:
                if legacy_pident >= 95.0 and legacy_cov >= 95.0:
                    legacy_hq += 1

            if updated_found and updated_pident is not None and updated_cov is not None:
                if updated_pident >= 95.0 and updated_cov >= 95.0:
                    updated_hq += 1

            if (
                legacy_found
                and updated_found
                and legacy_pident is not None
                and updated_pident is not None
                and legacy_cov is not None
                and updated_cov is not None
            ):
                identity_deltas.append(updated_pident - legacy_pident)
                coverage_deltas.append(updated_cov - legacy_cov)

    return MappingStats(
        total=total,
        legacy_hits=legacy_hits,
        updated_hits=updated_hits,
        both_hits=both_hits,
        legacy_only=legacy_only,
        updated_only=updated_only,
        no_hits=no_hits,
        legacy_hq=legacy_hq,
        updated_hq=updated_hq,
        identity_deltas=identity_deltas,
        coverage_deltas=coverage_deltas,
    )


def load_position_metrics(path: Path, limit: int = 5) -> List[PositionMetric]:
    metrics: List[PositionMetric] = []
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            coverage = int(float(row["coverage"]))
            mismatch_rate = float(row["mismatch_rate"])
            if coverage <= 0:
                continue
            metrics.append(
                PositionMetric(
                    position=int(row["position"]),
                    reference_base=row["reference_base"],
                    coverage=coverage,
                    mismatch_rate=mismatch_rate,
                )
            )
    metrics.sort(key=lambda item: item.mismatch_rate, reverse=True)
    return metrics[:limit]


def load_metric_summary(path: Path) -> dict[tuple[str, str], float]:
    summary: dict[tuple[str, str], float] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            summary[(row["dataset"], row["metric"])] = float(row["mean"])
    return summary


def load_variant_sites(path: Path) -> List[VariantSite]:
    sites: List[VariantSite] = []
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sites.append(
                VariantSite(
                    label=row["site"],
                    origin=row["origin"],
                    variant=row["variant"],
                    total_calls=int(float(row["total_calls"])),
                    variant_count=int(float(row["variant_count"])),
                    variant_percent=float(row["variant_percent"]),
                )
            )
    return sites


def format_delta(values: Iterable[float]) -> tuple[float, float, float]:
    data = list(values)
    if not data:
        return (0.0, 0.0, 0.0)
    return (mean(data), min(data), max(data))


def build_report(
    mapping: MappingStats,
    metric_summary: dict[tuple[str, str], float],
    top_positions: List[PositionMetric],
    variant_sites: List[VariantSite],
) -> str:
    mean_identity, min_identity, max_identity = format_delta(mapping.identity_deltas)
    mean_coverage, min_coverage, max_coverage = format_delta(mapping.coverage_deltas)

    lines: List[str] = []
    lines.append("Full operon mapping summary")
    lines.append("============================")
    lines.append(f"Assemblies scanned: {mapping.total}")
    lines.append(
        "Legacy hits (>= 90% identity, >= 80% coverage): "
        f"{mapping.legacy_hits}"
    )
    lines.append(
        "Updated hits (>= 90% identity, >= 80% coverage): "
        f"{mapping.updated_hits}"
    )
    lines.append(f"Hits in both references: {mapping.both_hits}")
    lines.append(f"Legacy-only hits: {mapping.legacy_only}")
    lines.append(f"Updated-only hits: {mapping.updated_only}")
    lines.append(f"Assemblies without hits: {mapping.no_hits}")
    lines.append("")
    lines.append("High-quality hits (>= 95% identity and coverage)")
    lines.append(f"  Legacy: {mapping.legacy_hq}")
    lines.append(f"  Updated: {mapping.updated_hq}")
    lines.append("")
    lines.append("Identity deltas (updated - legacy, hits in both references)")
    lines.append(f"  Mean: {mean_identity:.4f}")
    lines.append(f"  Min: {min_identity:.4f}")
    lines.append(f"  Max: {max_identity:.4f}")
    lines.append("")
    lines.append("Coverage deltas (updated - legacy, hits in both references)")
    lines.append(f"  Mean: {mean_coverage:.4f}")
    lines.append(f"  Min: {min_coverage:.4f}")
    lines.append(f"  Max: {max_coverage:.4f}")
    lines.append("")
    lines.append("Mean mismatch/deletion rates")
    lines.append(
        "  Legacy mismatch rate: "
        f"{metric_summary.get(('legacy', 'mismatch_rate'), 0.0):.6f}"
    )
    lines.append(
        "  Legacy deletion rate: "
        f"{metric_summary.get(('legacy', 'deletion_rate'), 0.0):.6f}"
    )
    lines.append(
        "  Updated mismatch rate: "
        f"{metric_summary.get(('updated', 'mismatch_rate'), 0.0):.6f}"
    )
    lines.append(
        "  Updated deletion rate: "
        f"{metric_summary.get(('updated', 'deletion_rate'), 0.0):.6f}"
    )
    lines.append("")
    lines.append("Top updated-reference mismatch positions (by mismatch rate)")
    for metric in top_positions:
        lines.append(
            f"  Pos {metric.position} (ref={metric.reference_base}): "
            f"coverage={metric.coverage} mismatch_rate={metric.mismatch_rate:.4f}"
        )
    lines.append("")

    if variant_sites:
        total_calls = variant_sites[0].total_calls
        variant_with_support = sum(1 for site in variant_sites if site.variant_count > 0)
        max_variant = max(site.variant_count for site in variant_sites)
        lines.append("Variant loci summarised in tables.md")
        lines.append(f"  Sites tracked: {len(variant_sites)}")
        lines.append(f"  Assemblies contributing allele calls per site: {total_calls}")
        lines.append(f"  Sites with non-zero variant support: {variant_with_support}")
        lines.append(f"  Maximum variant count across sites: {max_variant}")
        lines.append("  Variant allele counts:")
        for site in variant_sites:
            lines.append(
                f"    {site.label} ({site.origin}->{site.variant}): "
                f"{site.variant_count}/{site.total_calls} "
                f"({site.variant_percent:.4f}%)"
            )

    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mapping",
        type=Path,
        default=Path("output_full_run/full_operon_mapping.tsv"),
        help="Path to full_operon_mapping.tsv",
    )
    parser.add_argument(
        "--metric-summary",
        type=Path,
        default=Path("output_full_run/position_stats/position_metric_summary.tsv"),
        help="Path to position_metric_summary.tsv",
    )
    parser.add_argument(
        "--updated-position-stats",
        type=Path,
        default=Path("output_full_run/position_stats/updated_position_summary.tsv"),
        help="Path to updated_position_summary.tsv",
    )
    parser.add_argument(
        "--variant-summary",
        type=Path,
        default=Path("output_full_run/variant_site_summary.tsv"),
        help="Path to variant_site_summary.tsv",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output_full_run/manuscript_numbers.txt"),
        help="Destination for manuscript_numbers.txt",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=5,
        help="How many mismatch hot-spots to report",
    )
    args = parser.parse_args()

    mapping_stats = load_mapping_stats(args.mapping)
    metric_summary = load_metric_summary(args.metric_summary)
    top_positions = load_position_metrics(args.updated_position_stats, args.top)
    variant_sites = load_variant_sites(args.variant_summary)

    report = build_report(mapping_stats, metric_summary, top_positions, variant_sites)
    args.output.write_text(report)


if __name__ == "__main__":
    main()
