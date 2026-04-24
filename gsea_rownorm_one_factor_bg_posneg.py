#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import argparse
import numpy as np
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import textwrap


# =========================
# Basic utilities
# =========================

def read_auto(path):
    return pd.read_csv(path, sep=None, engine="python")


def detect_gene_col(df):
    for c in [
        "gene_symbol", "GeneSymbol", "SYMBOL",
        "gene", "Gene", "ENSEMBL", "ensembl", "Ensembl",
        "GeneID", "ID", "gene_ensembl"
    ]:
        if c in df.columns:
            return c
    return df.columns[0]


def load_orientation(map_path, factor):
    if map_path is None:
        return "+"
    ext = os.path.splitext(map_path)[-1].lower()
    if ext in [".yml", ".yaml"]:
        import yaml
        d = yaml.safe_load(open(map_path, "r"))
        o = str(d.get(factor, "+")).strip()
    else:
        m = pd.read_csv(map_path, sep=None, engine="python")
        fac_col = "factor" if "factor" in m.columns else m.columns[0]
        if "orientation" in m.columns:
            o = m.loc[m[fac_col].astype(str) == factor, "orientation"]
        else:
            o = m.loc[m[fac_col].astype(str) == factor, m.columns[1]]
        o = str(o.iloc[0]).strip() if len(o) > 0 else "+"
    return "+" if o in ["+", "plus", "pos", "1"] else "-"


def read_background_ids(background_path):
    bg = read_auto(background_path)
    gcol = detect_gene_col(bg)
    ids = bg[gcol].astype(str).dropna().drop_duplicates().tolist()
    return ids


def read_gmt_local(gmt_path):
    gene_sets = {}
    with open(gmt_path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            term = parts[0]
            genes = [g.strip() for g in parts[2:] if g.strip()]
            gene_sets[term] = list(dict.fromkeys(genes))
    if not gene_sets:
        raise RuntimeError(f"GMT を読めませんでした: {gmt_path}")
    return gene_sets


# =========================
# Library handling
# =========================

def ensure_gene_sets(library, universe_genes, targets=None):
    """
    library:
      - local GMT path
      - or gseapy library name
    """
    U = set(universe_genes)

    if os.path.isfile(library):
        enr = read_gmt_local(library)
    else:
        enr = gp.get_library(name=library, organism="Human")

    terms = {}
    for term, genes in enr.items():
        if targets:
            ok = any(t.upper() in str(term).upper() for t in targets)
            if not ok:
                continue
        inter = [g for g in genes if g in U]
        terms[str(term)] = list(dict.fromkeys(inter))

    if not terms:
        raise RuntimeError("フィルタ後の gene set が空です。library / targets / gene IDs を確認してください。")

    return terms


# =========================
# Row normalization
# =========================

def row_normalize(df_values: pd.DataFrame, mode: str, eps: float = 1e-12):
    X = df_values.to_numpy(dtype=float)

    if mode == "minmax":
        x_min = np.nanmin(X, axis=1, keepdims=True)
        x_max = np.nanmax(X, axis=1, keepdims=True)
        denom = np.maximum(x_max - x_min, eps)
        mm = (X - x_min) / denom
        mm = 2.0 * mm - 1.0
        return pd.DataFrame(mm, index=df_values.index, columns=df_values.columns)

    elif mode == "zscore":
        mu = np.nanmean(X, axis=1, keepdims=True)
        sd = np.nanstd(X, axis=1, keepdims=True)
        sd = np.maximum(sd, eps)
        z = (X - mu) / sd
        z = np.clip(z, -10, 10)
        return pd.DataFrame(z, index=df_values.index, columns=df_values.columns)

    elif mode == "rankpct":
        R = df_values.apply(lambda row: row.rank(method="average", ascending=False), axis=1)
        N = float(R.shape[1])
        pct = (N - R + 1.0) / N
        pct = 2.0 * pct - 1.0
        return pct

    elif mode == "none":
        return df_values.copy()

    else:
        raise ValueError("mode must be one of: minmax, zscore, rankpct, none")


def apply_ranking_mode(s: pd.Series, ranking_mode="raw"):
    s = s.astype(float)

    if ranking_mode == "raw":
        return s

    elif ranking_mode == "signed_rank":
        ranks = s.abs().rank(method="average", ascending=True)
        signed_ranks = ranks * np.sign(s)
        return signed_ranks

    else:
        raise ValueError("ranking_mode must be 'raw' or 'signed_rank'")


def make_stat_from_row_normalized(matrix_df, factor, orientation="+", mode="minmax",
                                  eps=1e-12, ranking_mode="raw"):
    df = matrix_df.copy()
    gcol = detect_gene_col(df)
    df.rename(columns={gcol: "gene_id"}, inplace=True)

    if factor not in df.columns:
        raise ValueError(f"列 '{factor}' が見つかりません。例: {list(df.columns[:10])}")

    values = df.drop(columns=["gene_id"])
    norm = row_normalize(values, mode=mode, eps=eps)

    s = norm[factor].astype(float)
    sgn = +1.0 if orientation in ["+", "plus", "pos", "1"] else -1.0
    s = s * sgn

    s = apply_ranking_mode(s, ranking_mode=ranking_mode)

    out = pd.DataFrame({
        "gene_id": df["gene_id"].astype(str).values,
        "stat": s.values
    })
    return out


# =========================
# Rank building
# =========================

def build_full_background_rank(stat_df, background_ids, jitter=1e-12, seed=42):
    """
    Fill missing genes with 0, then add tiny jitter only to exact zeros.
    Returns:
      ranked_series
      pre_jitter_series
    """
    s = pd.Series(0.0, index=pd.Index(background_ids, dtype=str), dtype=float)

    sub = (
        stat_df[["gene_id", "stat"]]
        .dropna(subset=["gene_id"])
        .drop_duplicates(subset=["gene_id"])
        .set_index("gene_id")["stat"]
        .astype(float)
    )

    sub_aligned = sub.reindex(s.index)
    mask = sub_aligned.notna()
    s.loc[mask] = sub_aligned[mask].values

    s_pre_jitter = s.copy()

    if jitter and jitter > 0:
        rng = np.random.default_rng(seed)
        zero_mask = (s.values == 0.0)
        if zero_mask.any():
            s.values[zero_mask] = rng.uniform(-jitter, jitter, size=zero_mask.sum())

    return s.sort_values(ascending=False), s_pre_jitter


# =========================
# Gene set filtering
# =========================

def filter_gene_sets_by_size_and_nonzero(gene_sets, nonzero_genes, min_size=10, max_size=5000, min_nonzero=0):
    keep = {}
    rows = []

    nonzero_genes = set(nonzero_genes)

    for term, genes in gene_sets.items():
        genes = list(dict.fromkeys(genes))
        n_universe = len(genes)
        n_nonzero = len(set(genes) & nonzero_genes)

        passed = (
            n_universe >= min_size and
            n_universe <= max_size and
            n_nonzero >= min_nonzero
        )

        rows.append({
            "Pathway": term,
            "n_genes_in_universe": n_universe,
            "n_nonzero_genes": n_nonzero,
            "keep": passed
        })

        if passed:
            keep[term] = genes

    qc = pd.DataFrame(rows).sort_values(
        ["keep", "n_nonzero_genes", "n_genes_in_universe"],
        ascending=[False, False, False]
    )

    return keep, qc


# =========================
# GSEA plot helpers
# =========================

def sanitize_filename(name, max_len=180):
    name = str(name)
    name = re.sub(r'[\\/:*?"<>|]+', "_", name)
    name = re.sub(r"\s+", "_", name).strip("_")
    if len(name) > max_len:
        name = name[:max_len]
    return name


def pick_sort_column(df):
    sort_candidates = ["FDR q-val", "fdr", "FDR", "qval", "padj", "NOM p-val"]
    for c in sort_candidates:
        if c in df.columns:
            return c
    return None


def pick_term_column(df):
    for c in ["Term", "Name", "Pathway"]:
        if c in df.columns:
            return c
    return df.columns[0]


def get_result_term_dict(pre, term):
    if term in pre.results:
        return pre.results[term]
    for k in pre.results.keys():
        if str(k) == str(term):
            return pre.results[k]
    raise KeyError(f"Term not found in pre.results: {term}")


def plot_gsea_simple(pre, term, outpath):
#def plot_gsea_simple(pre, term, outpath, ylim=(-1.0, 1.0)):
    res = get_result_term_dict(pre, term)

    es_profile = res.get("RES")
    hits = res.get("hits", [])
    nes = res.get("nes", res.get("NES", np.nan))

    if es_profile is None:
        print(f"[WARNING] No RES found for term: {term}")
        return

#    if ylim is not None:
#        ax.set_ylim(*ylim)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(es_profile, linewidth=2.5, color="#99cc33")
    ax.set_facecolor("white")
#    ax.spines["top"].set_visible(False)
#    ax.spines["right"].set_visible(False)
#    ax.set_ylim(-1.0, 1.0)

    if hits is not None and len(hits) > 0:
        ymin, ymax = ax.get_ylim()
        tick_bottom = ymin
        tick_top = ymin + (ymax - ymin) * 0.08
        for h in hits:
            ax.vlines(h, tick_bottom, tick_top, color="black", linewidth=0.6, alpha=0.6)

#    ax.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.axhline(0, color="gray", linewidth=0.8, alpha=0.6)

    wrapped_term = "\n".join(textwrap.wrap(str(term), width=50))
#    if np.isfinite(nes):
#        ax.set_title(f"{term}\nNES = {nes:.3f}")
#    else:
#        ax.set_title(f"{term}\nNES = NA")
    if nes is not None:
        ax.set_title(f"{wrapped_term}\nNES = {nes:.3f}", fontsize=11)
    else:
        ax.set_title(f"{wrapped_term}\nNES = NA", fontsize=11)

    ax.set_xlabel("Rank")
    ax.set_ylabel("Enrichment score")

#    fig.tight_layout()
#    fig.savefig(outpath, dpi=200)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)


def save_top_pos_neg_gsea_plots(pre, df_sorted, outdir, label, top_n=5):
    term_col = pick_term_column(df_sorted)

    if "NES" not in df_sorted.columns:
        print("[WARNING] NES column not found. Skip plotting.")
        return

    plot_dir = os.path.join(outdir, f"NES_only_plots_{label}")
    os.makedirs(plot_dir, exist_ok=True)

    df = df_sorted.copy()

    df_pos = df[df["NES"] > 0].sort_values("NES", ascending=False).head(top_n)
    df_neg = df[df["NES"] < 0].sort_values("NES", ascending=True).head(top_n)
    selected = pd.concat([df_pos, df_neg], axis=0)

    saved = []
    for i, (_, row) in enumerate(selected.iterrows(), start=1):
        term = row[term_col]
        nes = row["NES"]
        safe_term = sanitize_filename(term)
        sign = "POS" if nes > 0 else "NEG"
        outpath = os.path.join(plot_dir, f"{i:02d}_{sign}_{safe_term}_NES_only.png")

        try:
            plot_gsea_simple(pre, term, outpath)
            saved.append({
                "rank": i,
                "term": term,
                "NES": nes,
                "group": sign,
                "file": os.path.basename(outpath)
            })
        except Exception as e:
            print(f"[WARNING] Failed to plot {term}: {e}")

    if saved:
        pd.DataFrame(saved).to_csv(
            os.path.join(plot_dir, f"plot_manifest_pos_neg_{label}.csv"),
            index=False
        )


# =========================
# GSEA run
# =========================

def run_prerank_series(rnk_series, gene_sets, outdir, label,
                       min_size=10, max_size=5000, nperm=10000, threads=4,
                       weight=1.0, save_leading_edge=True, top_n_plots=5):
    os.makedirs(outdir, exist_ok=True)

    rnk = rnk_series.sort_values(ascending=False)

    pre = gp.prerank(
        rnk=rnk,
        gene_sets=gene_sets,
        outdir=os.path.join(outdir, f"report_{label}"),
        format="png",
        min_size=min_size,
        max_size=max_size,
        permutation_num=nperm,
        seed=42,
        no_plot=True,
        threads=threads,
        weight=weight
    )

    df = pre.res2d.copy()
    sort_col = pick_sort_column(df)
    df_sorted = df.sort_values(sort_col) if sort_col is not None else df.copy()

    gsea_csv = os.path.join(outdir, f"gsea_results_{label}.csv")
    df_sorted.to_csv(gsea_csv, index=False)

    rnk.to_csv(os.path.join(outdir, f"rnk_{label}.tsv"), sep="\t", header=False)

    if save_leading_edge:
        lead_cols = [c for c in df_sorted.columns if ("lead" in c.lower() or "ledge" in c.lower())]
        keep_cols = [c for c in ["Term", "Name", "NES", "NOM p-val", "FDR q-val"] if c in df_sorted.columns]
        keep_cols += lead_cols
        keep_cols = list(dict.fromkeys(keep_cols))
        if lead_cols:
            df_sorted[keep_cols].to_csv(
                os.path.join(outdir, f"gsea_leading_edge_{label}.csv"),
                index=False
            )

    save_top_pos_neg_gsea_plots(
        pre=pre,
        df_sorted=df_sorted,
        outdir=outdir,
        label=label,
        top_n=top_n_plots
    )

    return df_sorted


# =========================
# Main
# =========================

def main():
    ap = argparse.ArgumentParser(description="Completed improved rownorm GSEA for one factor across two groups")

    ap.add_argument("--factor", required=True)
    ap.add_argument("--matrix_A", required=True)
    ap.add_argument("--matrix_B", required=True)
    ap.add_argument("--label_A", default="K562")
    ap.add_argument("--label_B", default="NPC")

    ap.add_argument("--orient_map", default=None)
    ap.add_argument("--background", required=True)
    ap.add_argument("--library", required=True, help="local GMT path or gseapy library name")
    ap.add_argument("--targets", default="", help="comma-separated substring filter for pathway names")

    ap.add_argument("--min_size", type=int, default=10)
    ap.add_argument("--max_size", type=int, default=5000)
    ap.add_argument("--min_nonzero", type=int, default=0)

    ap.add_argument("--nperm", type=int, default=10000)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--weight", type=float, default=1.0)

    ap.add_argument("--rownorm", choices=["minmax", "zscore", "rankpct", "none"], default="minmax")
    ap.add_argument("--ranking_mode", choices=["raw", "signed_rank"], default="raw")
    ap.add_argument("--jitter", type=float, default=1e-12)
    ap.add_argument("--top_n_plots", type=int, default=5,
                    help="number of top positive and top negative NES pathways to save as NES-only plots")

    ap.add_argument("--save_leading_edge", action="store_true")
    ap.add_argument("--outdir", default="gsea_rownorm_one_factor_out")

    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    orientation = load_orientation(args.orient_map, args.factor)

    M_A = read_auto(args.matrix_A)
    M_B = read_auto(args.matrix_B)

    statA = make_stat_from_row_normalized(
        M_A, args.factor,
        orientation=orientation,
        mode=args.rownorm,
        ranking_mode=args.ranking_mode
    )
    statB = make_stat_from_row_normalized(
        M_B, args.factor,
        orientation=orientation,
        mode=args.rownorm,
        ranking_mode=args.ranking_mode
    )

    background_ids = read_background_ids(args.background)

    rnkA, preA = build_full_background_rank(statA, background_ids, jitter=args.jitter, seed=42)
    rnkB, preB = build_full_background_rank(statB, background_ids, jitter=args.jitter, seed=42)

    nonzeroA = set(preA[preA != 0].index)
    nonzeroB = set(preB[preB != 0].index)

    targets = [t.strip() for t in args.targets.split(",") if t.strip()] or None

    gene_sets_A_all = ensure_gene_sets(args.library, rnkA.index.tolist(), targets=targets)
    gene_sets_B_all = ensure_gene_sets(args.library, rnkB.index.tolist(), targets=targets)

    gene_sets_A, qcA = filter_gene_sets_by_size_and_nonzero(
        gene_sets_A_all, nonzeroA,
        min_size=args.min_size,
        max_size=args.max_size,
        min_nonzero=args.min_nonzero
    )
    gene_sets_B, qcB = filter_gene_sets_by_size_and_nonzero(
        gene_sets_B_all, nonzeroB,
        min_size=args.min_size,
        max_size=args.max_size,
        min_nonzero=args.min_nonzero
    )

    qcA.to_csv(os.path.join(args.outdir, f"gene_set_filter_qc_{args.label_A}.csv"), index=False)
    qcB.to_csv(os.path.join(args.outdir, f"gene_set_filter_qc_{args.label_B}.csv"), index=False)

    if len(gene_sets_A) == 0:
        raise RuntimeError(f"{args.label_A}: filtered gene sets are empty.")
    if len(gene_sets_B) == 0:
        raise RuntimeError(f"{args.label_B}: filtered gene sets are empty.")

    resA = run_prerank_series(
        rnkA, gene_sets_A, args.outdir, args.label_A,
        min_size=args.min_size,
        max_size=args.max_size,
        nperm=args.nperm,
        threads=args.threads,
        weight=args.weight,
        save_leading_edge=args.save_leading_edge,
        top_n_plots=args.top_n_plots
    )

    resB = run_prerank_series(
        rnkB, gene_sets_B, args.outdir, args.label_B,
        min_size=args.min_size,
        max_size=args.max_size,
        nperm=args.nperm,
        threads=args.threads,
        weight=args.weight,
        save_leading_edge=args.save_leading_edge,
        top_n_plots=args.top_n_plots
    )

    A = resA.copy()
    B = resB.copy()

    if "Term" in A.columns:
        A = A.rename(columns={"Term": "Pathway"})
    elif "Name" in A.columns:
        A = A.rename(columns={"Name": "Pathway"})
    else:
        A = A.rename(columns={A.columns[0]: "Pathway"})

    if "Term" in B.columns:
        B = B.rename(columns={"Term": "Pathway"})
    elif "Name" in B.columns:
        B = B.rename(columns={"Name": "Pathway"})
    else:
        B = B.rename(columns={B.columns[0]: "Pathway"})

    merged = pd.merge(
        A[["Pathway", "NES"]].rename(columns={"NES": f"NES_{args.label_A}"}),
        B[["Pathway", "NES"]].rename(columns={"NES": f"NES_{args.label_B}"}),
        on="Pathway",
        how="outer"
    )

    merged[f"deltaNES_{args.label_A}_minus_{args.label_B}"] = (
        merged[f"NES_{args.label_A}"] - merged[f"NES_{args.label_B}"]
    )

    merged.to_csv(
        os.path.join(args.outdir, f"deltaNES_{args.factor}_{args.label_A}_vs_{args.label_B}.csv"),
        index=False
    )

    print(
        f"[DONE] factor={args.factor} orientation={orientation} "
        f"rownorm={args.rownorm} ranking_mode={args.ranking_mode} "
        f"min_size={args.min_size} min_nonzero={args.min_nonzero} "
        f"top_n_plots={args.top_n_plots} "
        f"outdir={os.path.abspath(args.outdir)}"
    )


if __name__ == "__main__":
    main()
