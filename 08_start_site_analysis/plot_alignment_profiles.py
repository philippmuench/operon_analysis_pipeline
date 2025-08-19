#!/usr/bin/env python3
"""
Plot alignment profiles for upstream windows MSAs (e.g., ptsA_windows.aln.fa).

Outputs per-position profiles and a compact alignment heatmap:
- Top: conservation (1 - normalized Shannon entropy over A/C/G/T)
- Middle: gap frequency per position
- Bottom: stacked base composition (A/C/G/T) per position
- Separate PNG: alignment heatmap (subset of sequences) for visual inspection

Usage:
  python plot_alignment_profiles.py --msa output/ptsA_upstream/ptsA_windows.aln.fa \
    --outdir output/ptsA_upstream --max-rows 200
"""

import argparse
import os
from typing import Dict, List, Tuple, Optional

import numpy as np
import matplotlib.pyplot as plt
from Bio import AlignIO


BASES: Tuple[str, ...] = ("A", "C", "G", "T")
BASE_TO_IDX: Dict[str, int] = {b: i for i, b in enumerate(BASES)}


def compute_profiles(alignment) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute per-position base frequencies (A,C,G,T), gap frequency, and conservation.

    Returns:
      base_freqs: (L, 4) frequencies over A,C,G,T (ignoring gaps and N)
      gap_freq: (L,) fraction of sequences that are gaps at position
      conservation: (L,) 1 - H/Hmax, where H over A,C,G,T (gaps ignored), Hmax=log2(4)=2
    """
    L = alignment.get_alignment_length()
    N = len(alignment)

    base_counts = np.zeros((L, 4), dtype=float)
    gap_counts = np.zeros(L, dtype=float)

    for col in range(L):
        column = alignment[:, col]
        gc = 0
        for ch in column:
            ch = ch.upper()
            if ch == "-" or ch == ".":
                gap_counts[col] += 1
                continue
            if ch in BASE_TO_IDX:
                base_counts[col, BASE_TO_IDX[ch]] += 1
                gc += 1
            # ignore N/other

    # Base frequencies ignoring gaps
    with np.errstate(divide="ignore", invalid="ignore"):
        denom = base_counts.sum(axis=1, keepdims=True)
        base_freqs = np.divide(base_counts, denom, out=np.zeros_like(base_counts), where=denom > 0)

    # Gap frequency
    gap_freq = gap_counts / float(N)

    # Shannon entropy over A/C/G/T (ignore gaps) and normalized conservation
    eps = 1e-12
    H = -(base_freqs * np.log2(base_freqs + eps)).sum(axis=1)
    Hmax = np.log2(4.0)  # 2
    conservation = 1.0 - (H / Hmax)
    conservation = np.clip(conservation, 0.0, 1.0)

    return base_freqs, gap_freq, conservation


def plot_profiles(
    base_freqs: np.ndarray,
    gap_freq: np.ndarray,
    conservation: np.ndarray,
    out_png: str,
    title: str,
    save_pdf: bool = False,
    overlay_cols: Optional[Tuple[int, int, int]] = None,
    labels: Optional[Tuple[str, str]] = None,
) -> None:
    L = base_freqs.shape[0]
    x = np.arange(1, L + 1)

    fig, axes = plt.subplots(3, 1, figsize=(max(10, L / 12), 8), sharex=True, constrained_layout=True)

    # Top: conservation
    axes[0].plot(x, conservation, color="#1f77b4", lw=1.5)
    axes[0].set_ylabel("Conservation\n(1 - H/2)")
    axes[0].set_ylim(0, 1)
    axes[0].set_title(title)

    # Middle: gap frequency
    axes[1].plot(x, gap_freq, color="#d62728", lw=1.2)
    axes[1].set_ylabel("Gap frequency")
    axes[1].set_ylim(0, 1)

    # Bottom: stacked base composition
    colors = {"A": "#4daf4a", "C": "#377eb8", "G": "#ff7f00", "T": "#e41a1c"}
    bottom = np.zeros(L)
    for i, base in enumerate(BASES):
        axes[2].bar(x, base_freqs[:, i], bottom=bottom, width=1.0, color=colors[base], edgecolor="none", label=base)
        bottom += base_freqs[:, i]
    axes[2].set_ylabel("Base composition")
    axes[2].set_xlabel("Alignment position")
    axes[2].legend(loc="upper right", ncol=4, frameon=False)

    # Optional overlay of raw window region and start codon
    if overlay_cols is not None:
        up_col, start_col, end_col = overlay_cols
        for ax in axes:
            # Window span background
            ax.axvspan(max(1, up_col), min(L, end_col), color="#000000", alpha=0.08)
            
            # Region separators
            ax.axvline(max(1, start_col), color="#FF4444", lw=2.0, ls="-", alpha=0.8)
            ax.axvline(min(L, start_col + 3), color="#FF4444", lw=2.0, ls="-", alpha=0.8)
            
            # Start codon highlighting
            ax.axvspan(max(1, start_col), min(L, start_col + 3), color="#FFD54F", alpha=0.4)
        
        # Region labels on bottom axis
        upstream_center = (max(1, up_col) + max(1, start_col)) // 2
        downstream_center = (min(L, start_col + 3) + min(L, end_col)) // 2
        
        up_label = "UPSTREAM" if labels is None else labels[0]
        dn_label = "DOWNSTREAM" if labels is None else labels[1]
        axes[-1].text(upstream_center, 0.05, f"{up_label}", fontsize=8, color="#000000", 
                     ha="center", va="bottom", weight="bold",
                     bbox=dict(facecolor="lightblue", alpha=0.7, edgecolor="none", pad=1))
        axes[-1].text(max(1, start_col) + 1.5, 0.15, "START", fontsize=8, color="#000000", 
                     ha="center", va="bottom", weight="bold",
                     bbox=dict(facecolor="yellow", alpha=0.8, edgecolor="none", pad=1))
        axes[-1].text(downstream_center, 0.05, f"{dn_label}", fontsize=8, color="#000000", 
                     ha="center", va="bottom", weight="bold",
                     bbox=dict(facecolor="lightgreen", alpha=0.7, edgecolor="none", pad=1))

    fig.savefig(out_png, dpi=200)
    if save_pdf and out_png.endswith(".png"):
        fig.savefig(out_png[:-4] + ".pdf")
    plt.close(fig)


def plot_alignment_heatmap(alignment, out_png: str, max_rows: int = 200,
                           raw_windows: str = None, upstream_len: int = 120,
                           downstream_len: int = 60,
                           order_by: str = "consensus_distance", save_pdf: bool = False,
                           overlay_cols: Optional[Tuple[int, int, int]] = None,
                           full_gene: bool = False,
                           split_segments: bool = False) -> Optional[Tuple[int, int, int]]:
    L = alignment.get_alignment_length()
    N = len(alignment)
    rows = N if (max_rows is None or max_rows <= 0) else min(N, max_rows)

    # Sample evenly across sequences to get a representative subset
    idxs = np.linspace(0, N - 1, num=rows, dtype=int)
    # Cast numpy int to native int for Biopython indexer
    sub = [alignment[int(i)] for i in idxs]

    # Map A/C/G/T/- to ints for color mapping
    cmap_map = {"A": 0, "C": 1, "G": 2, "T": 3, "-": 4}
    arr = np.zeros((rows, L), dtype=int)
    for r, rec in enumerate(sub):
        s = str(rec.seq).upper()
        for c, ch in enumerate(s):
            if ch not in cmap_map:
                ch = "-" if ch in (".",) else "-"  # treat others as gap
            arr[r, c] = cmap_map[ch]

    # Optional ordering by consensus distance (default) before plotting
    order_idx = np.arange(rows)
    if order_by == "consensus_distance":
        cons = np.zeros(L, dtype=int)
        for c in range(L):
            counts = np.bincount(arr[:, c], minlength=5)
            cons[c] = int(np.argmax(counts))
        dists = (arr != cons).sum(axis=1)
        order_idx = np.argsort(dists, kind="mergesort")
        arr = arr[order_idx, :]
        sub = [sub[int(i)] for i in order_idx]

    # Define colors
    from matplotlib.colors import ListedColormap
    colors = ["#4daf4a", "#377eb8", "#ff7f00", "#e41a1c", "#cccccc"]  # A,C,G,T,-
    cmap = ListedColormap(colors)

    # Use a narrow side axis for start-codon class bar
    import matplotlib.gridspec as gridspec
    width_in = max(10, L / 12)
    height_in = min(40, max(4, rows / 250))
    fig = None
    ax = None
    ax_side = None

    # Defer drawing until after potential re-ordering based on start column
    def draw_matrix():
        nonlocal fig, ax, ax_side
        fig = plt.figure(figsize=(width_in, height_in), constrained_layout=True)
        import matplotlib.gridspec as gridspec  # local import to avoid top-level ref
        gs = gridspec.GridSpec(1, 2, width_ratios=[40, 1], wspace=0.05)
        ax = fig.add_subplot(gs[0, 0])
        ax_side = fig.add_subplot(gs[0, 1])
        im = ax.imshow(arr, aspect="auto", interpolation="nearest", cmap=cmap, origin="upper", zorder=0)
        ax.set_xlabel("Alignment position")
        ax.set_ylabel(f"Sequences (n={N}, shown={rows})")
        ax.set_title("Alignment heatmap (A/C/G/T/-)")
        return im

    # Create a legend proxy for base colors (A/C/G/T/Gaps)
    import matplotlib.patches as mpatches
    legend_patches = [
        mpatches.Patch(color=colors[0], label="A"),
        mpatches.Patch(color=colors[1], label="C"),
        mpatches.Patch(color=colors[2], label="G"),
        mpatches.Patch(color=colors[3], label="T"),
        mpatches.Patch(color=colors[4], label="Gap"),
    ]
    # We'll attach figure-level legends to avoid covering the heatmap

    # Optional: annotate start codon positions and codon types per row
    annotated = 0
    side_classes = np.zeros(rows, dtype=int)  # 0 unknown, 1 ATG, 2 GTG, 3 TTG
    med_up_col = None
    med_start_col = None
    med_end_col = None
    if raw_windows and os.path.exists(raw_windows):
        from Bio import SeqIO
        # Map record.id -> raw sequence
        raw_map = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(raw_windows, "fasta")}
        # Determine per-row aligned start column and codon
        codon_colors = {"ATG": "#2ca02c", "GTG": "#ff7f0e", "TTG": "#1f77b4"}
        codon_to_class = {"ATG": 1, "GTG": 2, "TTG": 3}
        start_cols_all: List[int] = []
        up_cols_all: List[int] = []
        end_cols_all: List[int] = []
        # Prepare overlay for subset rows
        ys: List[int] = []
        xs: List[int] = []
        cs: List[str] = []
        for r, rec in enumerate(sub):
            rid = rec.id
            raw = raw_map.get(rid)
            if not raw or len(raw) < upstream_len + 3:
                continue
            start_codon = raw[upstream_len: upstream_len + 3]
            # Map upstream_len (ungapped index) to aligned column index
            aligned_row = str(rec.seq).upper()
            ungapped = 0
            aligned_col = None
            up_col = None
            end_col = None
            for c in range(L):
                if aligned_row[c] != "-" and aligned_row[c] != ".":
                    if ungapped == 0 and up_col is None:
                        up_col = c
                    if ungapped == upstream_len and aligned_col is None:
                        aligned_col = c
                    if not full_gene:
                        if ungapped == upstream_len + 3 + downstream_len and end_col is None:
                            end_col = c
                    ungapped += 1
            if full_gene and end_col is None:
                # set to the last non-gap position in the aligned row
                for c in range(L - 1, -1, -1):
                    if aligned_row[c] != "-" and aligned_row[c] != ".":
                        end_col = c
                        break
            if aligned_col is None or up_col is None or end_col is None:
                continue
            start_cols_all.append(aligned_col)
            up_cols_all.append(up_col)
            end_cols_all.append(end_col)
            if start_codon in codon_colors:
                ys.append(r)
                xs.append(aligned_col)
                cs.append(codon_colors[start_codon])
                annotated += 1
                side_classes[r] = codon_to_class[start_codon]
        # Do not draw yet; we will draw after potential reordering
        if start_cols_all:
            med_col = int(np.median(start_cols_all))
        if up_cols_all and end_cols_all and start_cols_all:
            med_up_col = int(np.median(up_cols_all))
            med_end_col = int(np.median(end_cols_all))
            med_start_col = int(np.median(start_cols_all))
        # We previously added a scatter legend; now we rely on the side color strip + figure legend
        # Optional second ordering by start column (apply before drawing)
        if order_by == "start_column" and start_cols_all:
            aligned_cols_vec = np.full(rows, fill_value=10**9, dtype=int)
            for yi, xi in zip(ys, xs):
                aligned_cols_vec[yi] = xi
            order2 = np.argsort(aligned_cols_vec, kind="mergesort")
            inv_order = np.empty(rows, dtype=int)
            inv_order[order2] = np.arange(rows)
            arr = arr[order2, :]
            sub = [sub[int(i)] for i in order2]
            side_classes = side_classes[order2]
            # Remap scatter row indices to new order
            ys = [int(inv_order[y]) for y in ys]
        
        # Now draw matrix or split panels
        if split_segments and start_cols_all and up_cols_all and end_cols_all:
            # Build three-panel layout: Upstream | START | CDS, plus side bar
            upstream_w = max(1, med_start_col - med_up_col)
            start_w = 3
            cds_w = max(1, med_end_col - (med_start_col + 3))
            total_cols = upstream_w + start_w + cds_w
            # Figure and gridspec
            fig = plt.figure(figsize=(max(10, total_cols / 12), height_in), constrained_layout=True)
            import matplotlib.gridspec as gridspec
            gs = gridspec.GridSpec(1, 4, width_ratios=[upstream_w, start_w, cds_w, 1], wspace=0.05)
            ax_up = fig.add_subplot(gs[0, 0])
            ax_start = fig.add_subplot(gs[0, 1])
            ax_cds = fig.add_subplot(gs[0, 2])
            ax_side = fig.add_subplot(gs[0, 3])
            # Slices
            upstream_slice = arr[:, med_up_col:med_start_col]
            start_slice = arr[:, med_start_col:med_start_col + 3]
            cds_slice = arr[:, med_start_col + 3:med_end_col]
            ax_up.imshow(upstream_slice, aspect="auto", interpolation="nearest", cmap=cmap, origin="upper")
            ax_start.imshow(start_slice, aspect="auto", interpolation="nearest", cmap=cmap, origin="upper")
            ax_cds.imshow(cds_slice, aspect="auto", interpolation="nearest", cmap=cmap, origin="upper")
            # Labels/titles
            ax_up.set_title("Upstream")
            ax_start.set_title("START")
            ax_cds.set_title("CDS")
            ax_cds.set_xlabel("Alignment position (CDS region)")
            for a in (ax_up, ax_start, ax_cds):
                a.set_yticks([])
            ax_up.set_ylabel(f"Sequences (n={N}, shown={rows})")
            # Legend handled at figure level
            did_split = True
        else:
            im = draw_matrix()
            # Legend handled at figure level
            did_split = False
        
        # Scatter per-row start markers
        if not split_segments and xs and ys:
            ax.scatter(xs, ys, s=10, c=cs, marker="v", edgecolors="none", alpha=0.9, label="Start codon")
        # Draw a vertical line at median start column across subset
        if not split_segments and start_cols_all:
            ax.axvline(med_col, color="#000000", ls="--", lw=1.0, alpha=0.7)
            ax.text(med_col + 2, 3, "start", color="#000000", fontsize=8, va="bottom")
        # Draw overlays (window box and region labels)
        if not split_segments and med_up_col is not None and med_start_col is not None and med_end_col is not None:
            import matplotlib.patches as mpatches
            box = mpatches.Rectangle((med_up_col, -0.5), max(1, med_end_col - med_up_col), rows, fill=False,
                                     edgecolor="#111111", linewidth=2.5, linestyle="--", zorder=10)
            ax.add_patch(box)
            ax.axvspan(med_up_col, med_end_col, color="#000000", alpha=0.12, zorder=5)
            ax.axvline(med_start_col, color="#FF4444", lw=3.0, ls="-", zorder=15, alpha=0.8)
            ax.axvline(med_start_col + 3, color="#FF4444", lw=3.0, ls="-", zorder=15, alpha=0.8)
            upstream_center = (med_up_col + med_start_col) // 2
            downstream_center = (med_start_col + 3 + med_end_col) // 2
            up_text = "UPSTREAM\n(-%d bp)" % (upstream_len)
            dn_text = (
                "CDS\n(to gene end)" if full_gene else "DOWNSTREAM\n(+%d bp)" % (downstream_len)
            )
            ax.text(upstream_center, max(5, int(rows * 0.03)), up_text, 
                   fontsize=9, color="#000000", ha="center", va="bottom", zorder=16, weight="bold",
                   bbox=dict(facecolor="lightblue", alpha=0.7, edgecolor="none", pad=2))
            ax.text(med_start_col + 1.5, max(5, int(rows * 0.03)), "START\nCODON", 
                   fontsize=8, color="#000000", ha="center", va="bottom", zorder=16, weight="bold",
                   bbox=dict(facecolor="yellow", alpha=0.8, edgecolor="none", pad=1.5))
            ax.text(downstream_center, max(5, int(rows * 0.03)), dn_text, 
                   fontsize=9, color="#000000", ha="center", va="bottom", zorder=16, weight="bold",
                   bbox=dict(facecolor="lightgreen", alpha=0.7, edgecolor="none", pad=2))
            # Start codon band
            start_band = mpatches.Rectangle((med_start_col, -0.5), 3, rows, facecolor="#FFD54F",
                                            edgecolor="none", alpha=0.35, zorder=9)
            ax.add_patch(start_band)
    # Draw side column with start codon class
    from matplotlib.colors import ListedColormap
    side_cmap = ListedColormap(["#dddddd", "#2ca02c", "#ff7f0e", "#1f77b4"])  # unknown, ATG, GTG, TTG
    if ax_side is None:
        # Create minimal side axis if not already created (e.g., no raw_windows)
        fig = plt.figure(figsize=(width_in, height_in), constrained_layout=True)
        import matplotlib.gridspec as gridspec
        gs = gridspec.GridSpec(1, 2, width_ratios=[40, 1], wspace=0.05)
        ax = fig.add_subplot(gs[0, 0])
        ax_side = fig.add_subplot(gs[0, 1])
        ax.imshow(arr, aspect="auto", interpolation="nearest", cmap=cmap, origin="upper")
        ax.set_xlabel("Alignment position")
        ax.set_ylabel(f"Sequences (n={N}, shown={rows})")
        ax.set_title("Alignment heatmap (A/C/G/T/-)")
    ax_side.imshow(side_classes.reshape(rows, 1), aspect="auto", interpolation="nearest", cmap=side_cmap, origin="upper")
    ax_side.set_xticks([])
    ax_side.set_yticks([])
    ax_side.set_title("start codon type")
    # Print a small footer text with annotation count (attach to left-most axis available)
    # Attach footer text to a valid main axis
    if ax is not None:
        ax.text(5, rows - 5, f"annotated starts: {annotated}", fontsize=7, color="#333333", va="top")
    elif 'ax_up' in locals() and ax_up is not None:
        ax_up.text(2, rows - 5, f"annotated starts: {annotated}", fontsize=7, color="#333333", va="top")

    # Add figure-level legends: bases and start-codon types
    fig.legend(handles=legend_patches, loc="upper center", ncol=5, bbox_to_anchor=(0.5, 1.01), frameon=False)
    start_type_handles = [
        mpatches.Patch(color="#dddddd", label="Unknown"),
        mpatches.Patch(color="#2ca02c", label="ATG"),
        mpatches.Patch(color="#ff7f0e", label="GTG"),
        mpatches.Patch(color="#1f77b4", label="TTG"),
    ]
    fig.legend(handles=start_type_handles, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.02), frameon=False)

    # Use constrained_layout instead of tight_layout to avoid warnings with GridSpec
    fig.savefig(out_png, dpi=200)
    if save_pdf and out_png.endswith(".png"):
        fig.savefig(out_png[:-4] + ".pdf")
    plt.close(fig)

    # Return overlay columns for use in profiles plot
    if med_up_col is not None and med_start_col is not None and med_end_col is not None:
        return (med_up_col, med_start_col, med_end_col)
    # Or pass through provided overlay if given
    if overlay_cols is not None:
        return overlay_cols
    return None


def main() -> int:
    p = argparse.ArgumentParser(description="Plot alignment profiles for upstream windows MSA")
    p.add_argument("--msa", required=True, help="Path to aligned FASTA (e.g., *_windows.aln.fa)")
    p.add_argument("--outdir", default=None, help="Output directory (default: alongside MSA)")
    p.add_argument("--max-rows", type=int, default=200, help="Max sequences for heatmap subset (0/neg = all)")
    p.add_argument("--raw-windows", default=None, help="Raw windows FASTA (to annotate start codon)")
    p.add_argument("--upstream-len", type=int, default=120, help="Upstream length used when creating windows")
    p.add_argument("--downstream-len", type=int, default=60, help="Downstream length used when creating windows")
    p.add_argument("--order-by", choices=["consensus_distance", "start_column", "none"], default="consensus_distance",
                   help="Row ordering strategy")
    p.add_argument("--save-pdf", action="store_true", help="Also save PDF copies of plots")
    p.add_argument("--full-gene", dest="full_gene", action="store_true", help="If set, interpret windows as upstream_len upstream plus the entire CDS downstream of the start site.")
    p.add_argument("--split-segments", action="store_true", help="Split heatmap into Upstream | START | CDS panels with gaps between them.")
    args = p.parse_args()

    msa_path = args.msa
    if not os.path.exists(msa_path) or os.path.getsize(msa_path) == 0:
        print(f"Error: MSA not found or empty: {msa_path}")
        return 1

    outdir = args.outdir or os.path.dirname(msa_path) or "."
    os.makedirs(outdir, exist_ok=True)

    base = os.path.splitext(os.path.basename(msa_path))[0]
    title = f"{base}"

    print(f"Reading alignment: {msa_path}")
    alignment = AlignIO.read(msa_path, "fasta")
    print(f"  Sequences: {len(alignment)} | Length: {alignment.get_alignment_length()} bp")

    print("Computing profiles...")
    base_freqs, gap_freq, conservation = compute_profiles(alignment)

    out_profiles = os.path.join(outdir, f"{base}_profiles.png")
    print(f"Saving profiles plot: {out_profiles}")
    # Compute overlays via heatmap to align raw-window spans to aligned columns
    out_heatmap = os.path.join(outdir, f"{base}_heatmap.png")
    print(f"Saving heatmap: {out_heatmap}")
    overlay_cols = plot_alignment_heatmap(
        alignment,
        out_heatmap,
        max_rows=args.max_rows,
        raw_windows=args.raw_windows,
        upstream_len=args.upstream_len,
        downstream_len=args.downstream_len,
        order_by=args.order_by,
        save_pdf=args.save_pdf,
        full_gene=args.full_gene,
        split_segments=args.split_segments,
    )

    # Now plot profiles with the same overlay region
    labels = (f"UPSTREAM (-{args.upstream_len}bp)", f"DOWNSTREAM (+{args.downstream_len}bp)")
    if args.full_gene:
        labels = ("UPSTREAM", "CDS (to gene end)")
    plot_profiles(base_freqs, gap_freq, conservation, out_profiles, title, save_pdf=args.save_pdf, overlay_cols=overlay_cols, labels=labels)

    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


