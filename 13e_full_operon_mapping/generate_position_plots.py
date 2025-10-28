#!/usr/bin/env python3
"""Generate conservation/mismatch plots from combined position statistics."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional

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


def annotation_zoom_plot(
    df: pd.DataFrame,
    column: str,
    spans: Optional[List[Dict[str, float]]],
    target_label: str,
    title: str,
    ylabel: str,
    path: Path,
    color: str = "black",
    flank: int = 200,
    ylim: Optional[tuple[float, float]] = None,
):
    if not spans:
        return
    target = next((span for span in spans if span["label"].lower() == target_label.lower()), None)
    if not target:
        return

    min_clip = max(1, int(target["start"]) - flank)
    max_clip = min(int(df["position"].max()), int(target["end"]) + flank)
    subset = df[(df["position"] >= min_clip) & (df["position"] <= max_clip)].copy()
    if subset.empty:
        return

    display_label = target["label"].replace("_", " ").title()
    relevant_spans = [span for span in spans if span["end"] >= min_clip and span["start"] <= max_clip]

    line_plot(
        df=subset,
        y_cols=[column],
        ylabel=ylabel,
        title=title + f" (zoom on {display_label})",
        path=path,
        spans=relevant_spans,
        color=color,
        ylim=ylim,
    )
def codon_zoom_plot(
    df: pd.DataFrame,
    column: str,
    start: int,
    end: int,
    title: str,
    ylabel: str,
    path: Path,
    spans: Optional[List[Dict[str, float]]] = None,
    color: str = "black",
    ylim: Optional[tuple[float, float]] = None,
    pdf_path: Optional[Path] = None,
    figsize: tuple[float, float] = (12, 4),
    pdf_figsize: Optional[tuple[float, float]] = None,
    plot_kind: str = "line",
    legend_label: Optional[str] = None,
):
    subset = df[(df["position"] >= start) & (df["position"] <= end)].copy()
    if subset.empty:
        return

    fig, ax = plt.subplots(1, 1, figsize=figsize)

    values = subset[column]
    label = legend_label or ylabel
    if plot_kind == "bar":
        ax.bar(subset["position"], values, width=0.8, color=color, align="center", label=label)
    else:
        ax.plot(subset["position"], values, color=color, linewidth=1, label=label)

    ax.set_ylabel(ylabel)
    ax.set_title(title + f" (positions {start}-{end})")
    if ylim is not None:
        ax.set_ylim(*ylim)
    elif plot_kind == "bar":
        y_max = values.max() if not values.empty else 0
        ax.set_ylim(0, max(1, y_max * 1.1))

    if spans:
        y_min, y_max = ax.get_ylim()
        y_text = y_max - (y_max - y_min) * 0.05
        for span in spans:
            if span["end"] < start or span["start"] > end:
                continue
            span_start = max(span["start"], start) - 0.5
            span_end = min(span["end"], end) + 0.5
            ax.axvspan(span_start, span_end, color=span["color"], alpha=0.15)
            ax.text(
                (span_start + span_end) / 2,
                y_text,
                span["label"],
                rotation=90,
                ha="center",
                va="top",
                fontsize=8,
            )

    ax.set_xlim(start - 0.5, end + 0.5)

    positions = subset["position"].astype(int).tolist()
    bases = subset["reference_base"].astype(str).str.upper().tolist()
    ax.set_xticks(positions)
    ax.set_xticklabels(bases, fontsize=5)
    ax.tick_params(axis="x", rotation=90, labelsize=5, pad=2)
    ax.set_xlabel("Reference base")

    top_ax = ax.twiny()
    top_ax.set_xlim(ax.get_xlim())
    tick_step = max(1, (end - start) // 10)
    top_ticks = list(range(start, end + 1, tick_step))
    if top_ticks[-1] != end:
        top_ticks.append(end)
    top_ax.set_xticks(top_ticks)
    top_ax.set_xticklabels([str(tick) for tick in top_ticks], fontsize=7)
    top_ax.set_xlabel("Operon position (nt)")

    fig.tight_layout()
    if legend_label:
        ax.legend(loc="upper right", fontsize=8)
    fig.savefig(path)

    if pdf_path:
        original_size = fig.get_size_inches()
        new_size = pdf_figsize or original_size
        fig.set_size_inches(new_size[0], new_size[1])
        fig.tight_layout()
        fig.savefig(pdf_path)
        fig.set_size_inches(*original_size)

    plt.close(fig)


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
    if "insertion_rate_legacy" in df:
        df["insertion_rate_legacy_pct"] = df["insertion_rate_legacy"] * 100
    else:
        df["insertion_rate_legacy_pct"] = 0.0
    if "insertion_rate_updated" in df:
        df["insertion_rate_updated_pct"] = df["insertion_rate_updated"] * 100
    else:
        df["insertion_rate_updated_pct"] = 0.0

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
    if "insertion_rate_legacy_pct" in df:
        line_plot(
            df=df,
            y_cols=["insertion_rate_legacy_pct"],
            ylabel="Insertion rate (%)",
            title=f"Insertion rate along operon ({args.legacy_label})",
            path=args.output_dir / "legacy_insertion_rate.png",
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
    annotation_zoom_plot(
        df=df,
        column="deletion_rate_legacy_pct",
        spans=legacy_spans,
        target_label="pribnow_box",
        title=f"Deletion rate along operon ({args.legacy_label})",
        ylabel="Deletion (gap) rate (%)",
        path=args.output_dir / "legacy_deletion_rate_pribnow_zoom.png",
        ylim=(0, 100),
    )
    codon_zoom_plot(
        df=df,
        column="deletion_count_legacy",
        start=7500,
        end=7600,
        title=f"Deletion rate along operon ({args.legacy_label})",
        ylabel="Deletion gaps (# assemblies)",
        path=args.output_dir / "legacy_deletion_rate_pribnow_zoom_detail.png",
        spans=legacy_spans,
        pdf_path=args.output_dir / "legacy_deletion_rate_pribnow_zoom_detail.pdf",
        pdf_figsize=(18, 3.5),
        plot_kind="bar",
        legend_label="Deletion gaps",
    )
    if "insertion_count_legacy" in df:
        codon_zoom_plot(
            df=df,
            column="insertion_count_legacy",
            start=7500,
            end=7600,
            title=f"Insertion rate along operon ({args.legacy_label})",
            ylabel="Insertion events (# assemblies)",
            path=args.output_dir / "legacy_insertion_rate_pribnow_zoom_detail.png",
            spans=legacy_spans,
            pdf_path=args.output_dir / "legacy_insertion_rate_pribnow_zoom_detail.pdf",
            pdf_figsize=(18, 3.5),
            plot_kind="bar",
            legend_label="Insertion events",
        )
    codon_zoom_plot(
        df=df,
        column="mismatch_rate_legacy_pct",
        start=7500,
        end=7600,
        title=f"Mismatch rate along operon ({args.legacy_label})",
        ylabel="Mismatch rate (%)",
        path=args.output_dir / "legacy_mismatch_rate_pribnow_zoom_detail.png",
        spans=legacy_spans,
        pdf_path=args.output_dir / "legacy_mismatch_rate_pribnow_zoom_detail.pdf",
        pdf_figsize=(18, 3.5),
        plot_kind="bar",
        legend_label="Mismatch rate",
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
    if "insertion_rate_updated_pct" in df:
        line_plot(
            df=df,
            y_cols=["insertion_rate_updated_pct"],
            ylabel="Insertion rate (%)",
            title=f"Insertion rate along operon ({args.updated_label})",
            path=args.output_dir / "updated_insertion_rate.png",
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
    annotation_zoom_plot(
        df=df,
        column="deletion_rate_updated_pct",
        spans=updated_spans,
        target_label="pribnow_box",
        title=f"Deletion rate along operon ({args.updated_label})",
        ylabel="Deletion (gap) rate (%)",
        path=args.output_dir / "updated_deletion_rate_pribnow_zoom.png",
        ylim=(0, 100),
    )
    codon_zoom_plot(
        df=df,
        column="deletion_count_updated",
        start=7500,
        end=7600,
        title=f"Deletion rate along operon ({args.updated_label})",
        ylabel="Deletion gaps (# assemblies)",
        path=args.output_dir / "updated_deletion_rate_pribnow_zoom_detail.png",
        spans=updated_spans,
        pdf_path=args.output_dir / "updated_deletion_rate_pribnow_zoom_detail.pdf",
        pdf_figsize=(18, 3.5),
        plot_kind="bar",
        legend_label="Deletion gaps",
    )
    if "insertion_count_updated" in df:
        codon_zoom_plot(
            df=df,
            column="insertion_count_updated",
            start=7500,
            end=7600,
            title=f"Insertion rate along operon ({args.updated_label})",
            ylabel="Insertion events (# assemblies)",
            path=args.output_dir / "updated_insertion_rate_pribnow_zoom_detail.png",
            spans=updated_spans,
            pdf_path=args.output_dir / "updated_insertion_rate_pribnow_zoom_detail.pdf",
            pdf_figsize=(18, 3.5),
            plot_kind="bar",
            legend_label="Insertion events",
        )
    codon_zoom_plot(
        df=df,
        column="mismatch_rate_updated_pct",
        start=7500,
        end=7600,
        title=f"Mismatch rate along operon ({args.updated_label})",
        ylabel="Mismatch rate (%)",
        path=args.output_dir / "updated_mismatch_rate_pribnow_zoom_detail.png",
        spans=updated_spans,
        pdf_path=args.output_dir / "updated_mismatch_rate_pribnow_zoom_detail.pdf",
        pdf_figsize=(18, 3.5),
        plot_kind="bar",
        legend_label="Mismatch rate",
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
