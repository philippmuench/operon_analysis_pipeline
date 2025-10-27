#!/usr/bin/env python3
"""Generate conservation/mismatch plots from combined position statistics."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ColorPalette = [
    "#4C72B0",
    "#DD8452",
    "#55A868",
    "#C44E52",
    "#8172B2",
    "#937860",
    "#DA8BC3",
    "#8C8C8C",
    "#CCB974",
    "#64B5CD",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("output_full_run/position_stats/combined_position_summary.tsv"),
        help="Path to combined position summary TSV",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output_full_run/position_stats/plots"),
        help="Directory for generated plots",
    )
    parser.add_argument(
        "--legacy-label",
        default="Legacy",
        help="Label for legacy dataset",
    )
    parser.add_argument(
        "--updated-label",
        default="Updated",
        help="Label for updated dataset",
    )
    parser.add_argument(
        "--legacy-annotation",
        type=Path,
        help="Optional TSV with gene/promoter spans for the legacy reference",
    )
    parser.add_argument(
        "--updated-annotation",
        type=Path,
        help="Optional TSV with gene/promoter spans for the updated reference",
    )
    parser.add_argument(
        "--deletion-threshold",
        type=float,
        default=1.5,
        help="Threshold (in %) for generating zoomed deletion-rate plots",
    )
    return parser.parse_args()


def line_plot(
    df: pd.DataFrame,
    y_cols: list[str],
    ylabel: str,
    title: str,
    path: Path,
    spans: Optional[List[Dict[str, float]]] = None,
    color: str = "black",
    ylim: Optional[tuple[float, float]] = None,
):
    plt.figure(figsize=(10, 4))
    for col in y_cols:
        if col in df:
            plt.plot(df["position"], df[col], color=color, linewidth=1)
    plt.xlabel("Position along operon (nt)")
    plt.ylabel(ylabel)
    plt.title(title)
    if ylim is not None:
        plt.ylim(*ylim)

    ax = plt.gca()
    if not df.empty:
        x_min = float(df["position"].min())
        x_max = float(df["position"].max())
        ax.set_xlim(x_min, x_max)

    if spans:
        y_min, y_max = ax.get_ylim()
        y_text = y_max - (y_max - y_min) * 0.05
        for span in spans:
            ax.axvspan(span["start"], span["end"], color=span["color"], alpha=0.15)
            ax.text(
                (span["start"] + span["end"]) / 2,
                y_text,
                span["label"],
                rotation=90,
                ha="center",
                va="top",
                fontsize=8,
            )

    plt.tight_layout()
    plt.savefig(path)
    plt.close()


def zoom_plot(
    df: pd.DataFrame,
    column: str,
    threshold_pct: float,
    title: str,
    ylabel: str,
    path: Path,
    spans: Optional[List[Dict[str, float]]] = None,
    color: str = "black",
    margin: int = 50,
    ylim: Optional[tuple[float, float]] = None,
):
    mask = df[column] > threshold_pct
    if not mask.any():
        return
    min_pos = int(df.loc[mask, "position"].min())
    max_pos = int(df.loc[mask, "position"].max())
    min_clip = max(1, min_pos - margin)
    max_clip = min(int(df["position"].max()), max_pos + margin)
    subset = df[(df["position"] >= min_clip) & (df["position"] <= max_clip)].copy()
    if subset.empty:
        return
    line_plot(
        df=subset,
        y_cols=[column],
        ylabel=ylabel,
        title=title + f" (positions {min_clip}-{max_clip})",
        path=path,
        spans=[span for span in (spans or []) if span["end"] >= min_clip and span["start"] <= max_clip],
        color=color,
        ylim=ylim,
    )


def load_annotations(path: Optional[Path]) -> List[Dict[str, float]]:
    if not path or not path.exists():
        return []
    df = pd.read_csv(path, sep="\t")
    spans: List[Dict[str, float]] = []
    colors = ColorPalette
    for idx, row in df.iterrows():
        label = row.get("gene") or row.get("locus_tag") or row.get("label") or f"feature_{idx}"
        lower_label = label.lower() if label else ""
        if lower_label in {"operon_promoter", "promoter"}:
            continue
        start = float(row["start"])
        end = float(row["end"]) or start
        spans.append(
            {
                "label": label,
                "start": start,
                "end": end,
                "color": colors[idx % len(colors)],
            }
        )
    return spans


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.input, sep="\t")
    if df.empty:
        raise SystemExit(f"No data found in {args.input}")

    df.sort_values("position", inplace=True)
    df.fillna(0, inplace=True)

    total_legacy = df["coverage_legacy"].max() or 1
    total_updated = df["coverage_updated"].max() or 1
    df["coverage_legacy_pct"] = df["coverage_legacy"] / total_legacy * 100
    df["coverage_updated_pct"] = df["coverage_updated"] / total_updated * 100
    df["mismatch_rate_legacy_pct"] = df["mismatch_rate_legacy"] * 100
    df["mismatch_rate_updated_pct"] = df["mismatch_rate_updated"] * 100
    df["deletion_rate_legacy_pct"] = df["deletion_rate_legacy"] * 100
    df["deletion_rate_updated_pct"] = df["deletion_rate_updated"] * 100

    args.output_dir.mkdir(parents=True, exist_ok=True)

    legacy_spans = load_annotations(args.legacy_annotation)
    updated_spans = load_annotations(args.updated_annotation)

    # Legacy plots
    line_plot(
        df=df,
        y_cols=["mismatch_rate_legacy_pct"],
        ylabel="Mismatch rate (%)",
        title=f"Mismatch rate along operon ({args.legacy_label})",
        path=args.output_dir / "legacy_mismatch_rate.png",
        spans=legacy_spans,
        ylim=(0, 100),
    )

    line_plot(
        df=df,
        y_cols=["deletion_rate_legacy_pct"],
        ylabel="Deletion (gap) rate (%)",
        title=f"Deletion rate along operon ({args.legacy_label})",
        path=args.output_dir / "legacy_deletion_rate.png",
        spans=legacy_spans,
        ylim=(0, 100),
    )
    zoom_plot(
        df=df,
        column="deletion_rate_legacy_pct",
        threshold_pct=args.deletion_threshold,
        title=f"Deletion rate along operon ({args.legacy_label})",
        ylabel="Deletion (gap) rate (%)",
        path=args.output_dir / "legacy_deletion_rate_zoom.png",
        spans=legacy_spans,
        ylim=(0, 100),
    )

    line_plot(
        df=df,
        y_cols=["coverage_legacy_pct"],
        ylabel="Coverage (% assemblies)",
        title=f"Coverage along operon ({args.legacy_label})",
        path=args.output_dir / "legacy_coverage.png",
        spans=legacy_spans,
        ylim=(0, 100),
    )

    # Updated plots
    line_plot(
        df=df,
        y_cols=["mismatch_rate_updated_pct"],
        ylabel="Mismatch rate (%)",
        title=f"Mismatch rate along operon ({args.updated_label})",
        path=args.output_dir / "updated_mismatch_rate.png",
        spans=updated_spans,
        ylim=(0, 100),
    )

    line_plot(
        df=df,
        y_cols=["deletion_rate_updated_pct"],
        ylabel="Deletion (gap) rate (%)",
        title=f"Deletion rate along operon ({args.updated_label})",
        path=args.output_dir / "updated_deletion_rate.png",
        spans=updated_spans,
        ylim=(0, 100),
    )
    zoom_plot(
        df=df,
        column="deletion_rate_updated_pct",
        threshold_pct=args.deletion_threshold,
        title=f"Deletion rate along operon ({args.updated_label})",
        ylabel="Deletion (gap) rate (%)",
        path=args.output_dir / "updated_deletion_rate_zoom.png",
        spans=updated_spans,
        ylim=(0, 100),
    )

    line_plot(
        df=df,
        y_cols=["coverage_updated_pct"],
        ylabel="Coverage (% assemblies)",
        title=f"Coverage along operon ({args.updated_label})",
        path=args.output_dir / "updated_coverage.png",
        spans=updated_spans,
        ylim=(0, 100),
    )


if __name__ == "__main__":
    main()
