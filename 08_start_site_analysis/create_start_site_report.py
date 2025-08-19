#!/usr/bin/env python3
"""
Create an interpretable report (plots + HTML) from start_site_summary.tsv.
"""

import os
import argparse
import pandas as pd
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import seaborn as sns  # optional
    sns.set_context("talk")
    sns.set_style("whitegrid")
except Exception:
    sns = None


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def load_data(tsv_path: str) -> pd.DataFrame:
    df = pd.read_csv(tsv_path, sep="\t")
    # Coerce numeric fields
    numeric_cols = [
        "start","end","upstream_offset_nt","upstream_rbs_spacer_nt","upstream_rbs_mismatches"
    ]
    for c in numeric_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    # Booleans
    for c in ["upstream_candidate_found","upstream_has_rbs"]:
        if c in df.columns:
            df[c] = df[c].astype(bool)
    # Categorical ordering
    if "classification" in df.columns:
        order = ["Alt_upstream_with_RBS","Upstream_no_RBS","No_upstream","Partial_gene"]
        df["classification"] = pd.Categorical(df["classification"], categories=order, ordered=True)
    return df


def plot_classification_by_gene(df: pd.DataFrame, out_png: str) -> None:
    counts = df.groupby(["gene","classification"]).size().reset_index(name="n")
    pivot = counts.pivot(index="gene", columns="classification", values="n").fillna(0)
    pivot = pivot.loc[pivot.sum(axis=1).sort_values(ascending=False).index]
    ax = pivot.plot(kind="bar", stacked=True, figsize=(14,6), colormap="tab20")
    ax.set_ylabel("Genomes")
    ax.set_xlabel("Gene")
    ax.set_title("Start-site classification per gene")
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()


def plot_upstream_offset_distribution(df: pd.DataFrame, out_png: str) -> None:
    d = df[df["upstream_offset_nt"].notna() & (df["upstream_offset_nt"] > 0)]
    if d.empty:
        return
    plt.figure(figsize=(12,5))
    if sns is not None:
        sns.histplot(d, x="upstream_offset_nt", bins=40, hue="gene", multiple="stack")
    else:
        for g, gdf in d.groupby("gene"):
            plt.hist(gdf["upstream_offset_nt"], bins=40, alpha=0.5, label=g)
        plt.legend(fontsize=8)
    plt.xlabel("Upstream offset (nt)")
    plt.ylabel("Count")
    plt.title("Distribution of upstream start offsets")
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()


def plot_rbs_rate_by_gene(df: pd.DataFrame, out_png: str) -> None:
    # Among cases with an upstream candidate, share that have RBS detected
    cand = df[df["upstream_candidate_found"]]
    if cand.empty:
        return
    rate = cand.groupby("gene")["upstream_has_rbs"].mean().reset_index(name="rbs_rate")
    rate = rate.sort_values("rbs_rate", ascending=False)
    plt.figure(figsize=(12,5))
    plt.bar(rate["gene"], rate["rbs_rate"], color="#4C72B0")
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("RBS rate among upstream candidates")
    plt.ylim(0,1)
    plt.title("RBS presence rate by gene")
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()


def plot_start_codon_distribution(df: pd.DataFrame, out_png: str) -> None:
    # Distribution of called start codons by gene
    counts = df.groupby(["gene","start_codon"]).size().reset_index(name="n")
    pivot = counts.pivot(index="gene", columns="start_codon", values="n").fillna(0)
    ax = pivot.plot(kind="bar", stacked=True, figsize=(14,6))
    ax.set_ylabel("Genomes")
    ax.set_xlabel("Gene")
    ax.set_title("Called start codon distribution by gene")
    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    plt.close()


def write_top_tables(df: pd.DataFrame, out_dir: str) -> dict:
    stats = {}
    # Alt_upstream_with_RBS fraction per gene
    total = df.groupby("gene").size().rename("total")
    alt_rbs = df[df["classification"] == "Alt_upstream_with_RBS"].groupby("gene").size().rename("alt_rbs")
    tab = pd.concat([total, alt_rbs], axis=1).fillna(0)
    tab["frac_alt_rbs"] = tab["alt_rbs"] / tab["total"].replace(0, np.nan)
    tab = tab.sort_values("frac_alt_rbs", ascending=False)
    path = os.path.join(out_dir, "top_genes_alt_upstream_with_rbs.tsv")
    tab.to_csv(path, sep="\t")
    stats["top_genes_table"] = path
    return stats


def write_html(report_dir: str, fig_paths: dict, table_paths: dict) -> str:
    html_path = os.path.join(report_dir, "report.html")
    with open(html_path, "w") as f:
        f.write("<html><head><meta charset='utf-8'><title>Start-site analysis report</title></head><body>")
        f.write("<h2>Start-site analysis report</h2>")
        for title, img in [
            ("Classification per gene", fig_paths.get("class_per_gene")),
            ("Upstream offset distribution", fig_paths.get("offset_dist")),
            ("RBS rate by gene", fig_paths.get("rbs_rate")),
            ("Called start codon distribution", fig_paths.get("start_codon_dist")),
        ]:
            if img and os.path.exists(img):
                f.write(f"<h3>{title}</h3><img src='{os.path.basename(img)}' style='max-width:1200px;'><br/>")
        if table_paths.get("top_genes_table"):
            rel = os.path.basename(table_paths["top_genes_table"])
            f.write(f"<h3>Top genes with Alt_upstream_with_RBS</h3><p><a href='{rel}'>{rel}</a></p>")
        f.write("</body></html>")
    return html_path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input_tsv", default=os.path.join("output", "start_site_summary.tsv"))
    ap.add_argument("--out_dir", default=os.path.join("output", "report"))
    args = ap.parse_args()

    ensure_dir(args.out_dir)

    df = load_data(args.input_tsv)
    fig_paths = {}

    p1 = os.path.join(args.out_dir, "classification_per_gene.png")
    plot_classification_by_gene(df, p1)
    fig_paths["class_per_gene"] = p1

    p2 = os.path.join(args.out_dir, "upstream_offset_distribution.png")
    plot_upstream_offset_distribution(df, p2)
    fig_paths["offset_dist"] = p2

    p3 = os.path.join(args.out_dir, "rbs_rate_by_gene.png")
    plot_rbs_rate_by_gene(df, p3)
    fig_paths["rbs_rate"] = p3

    p4 = os.path.join(args.out_dir, "start_codon_distribution.png")
    plot_start_codon_distribution(df, p4)
    fig_paths["start_codon_dist"] = p4

    tables = write_top_tables(df, args.out_dir)
    html_path = write_html(args.out_dir, fig_paths, tables)
    print(f"Wrote report: {html_path}")


if __name__ == "__main__":
    main()


