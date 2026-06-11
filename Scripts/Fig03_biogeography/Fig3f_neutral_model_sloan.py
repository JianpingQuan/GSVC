# GSVC_github publication script
# =============================================================================
# GSVC paper figure script (standalone copy for publication reproducibility)
# Regenerate: python paper_figure_scripts/generate_paper_scripts.py
# =============================================================================
from __future__ import annotations
import os
import sys
from pathlib import Path

PHAGE_ROOT = Path(os.environ.get("PHAGE_PROJECT_ROOT", Path(__file__).resolve().parents[2])).expanduser().resolve()
if not (PHAGE_ROOT / "FINAL_all_projects_Count_matrix.tsv").exists():
    PHAGE_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = PHAGE_ROOT / "Scripts" / "output"
OUT_DIR.mkdir(parents=True, exist_ok=True)
os.chdir(PHAGE_ROOT)


#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def norm_missing(x) -> str:
    if pd.isna(x):
        return "NA"
    s = str(x).strip()
    if s == "" or s.lower() in {"na", "nan", "none"}:
        return "NA"
    return s


def set_pub_rcparams():
    import matplotlib as mpl

    mpl.rcParams.update(
        {
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.family": "sans-serif",
            "font.sans-serif": [
                "Arial",
                "Helvetica",
                "Arial Unicode MS",
                "DejaVu Sans",
                "Liberation Sans",
            ],
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 11,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "figure.dpi": 300,
            "savefig.dpi": 300,
        }
    )


def fig_permanova_schemeA(outdir: Path, schemeA_combined: Path):
    import matplotlib.pyplot as plt
    import seaborn as sns

    d = pd.read_csv(schemeA_combined, sep="\t")
    # Normalize column names across possible outputs
    if "Pr(>F)" in d.columns:
        d = d.rename(columns={"Pr(>F)": "pvalue"})
    if "R2" not in d.columns:
        raise ValueError("schemeA combined must have R2 column")

    # Keep only key terms
    keep_terms = {"Country", "Gut_location", "Treatment", "Project"}
    d = d[d["term"].isin(keep_terms)].copy()
    d["model"] = d["model"].astype(str)

    # Order: without_project shows Country; with_project shows Project
    term_order = ["Country", "Gut_location", "Treatment", "Project"]
    d["term"] = pd.Categorical(d["term"], categories=term_order, ordered=True)

    # Barplot
    fig, ax = plt.subplots(1, 1, figsize=(8, 3.8), constrained_layout=True)
    sns.barplot(data=d, x="term", y="R2", hue="model", ax=ax, palette={"without_project": "#4C72B0", "with_project": "#55A868"})
    ax.set_title("PERMANOVA (adonis2, by=margin) on CLR-PC space (K=20)")
    ax.set_xlabel("")
    ax.set_ylabel("Marginal R²")
    ax.legend(title="", frameon=False, loc="upper right")
    ax.set_ylim(0, max(0.03, float(d["R2"].max()) * 1.25))
    for lab in ax.get_xticklabels():
        lab.set_rotation(30)
        lab.set_horizontalalignment("right")

    for ext in ["pdf", "png"]:
        fig.savefig(outdir / f"Fig5_permanova_schemeA_marginalR2.{ext}", bbox_inches="tight")
    plt.close(fig)


def neutral_sloan_fit(mean_abd: np.ndarray, occ: np.ndarray, n_reads: float) -> tuple[float, float, np.ndarray]:
    """
    Very lightweight Sloan-like fit:
    Use a logistic curve occ ~ 1 - exp(-m * mean_abd * n_reads)
    Fit m by grid search on log10(m).
    Returns best_m, rmse, pred_occ
    """
    # keep valid
    mgrid = np.logspace(-6, 2, 200)
    best_m = mgrid[0]
    best_rmse = 1e9
    best_pred = None
    for m in mgrid:
        pred = 1.0 - np.exp(-m * mean_abd * n_reads)
        rmse = float(np.sqrt(np.mean((pred - occ) ** 2)))
        if rmse < best_rmse:
            best_rmse = rmse
            best_m = float(m)
            best_pred = pred
    return best_m, best_rmse, best_pred  # type: ignore


def fig_neutral_model_top30000(outdir: Path, counts_top: Path):
    import matplotlib.pyplot as plt
    import seaborn as sns

    # pvca_Y_count.tsv: samples x features, first col sample_id
    df = pd.read_csv(counts_top, sep="\t")
    X = df.iloc[:, 1:].to_numpy(dtype=np.float64)
    # relative abundance
    lib = X.sum(axis=1, keepdims=True)
    lib[lib == 0] = 1.0
    rel = X / lib
    mean_abd = rel.mean(axis=0)
    occ = (X > 0).mean(axis=0)
    # effective reads scale (median library)
    n_reads = float(np.median(lib))

    # remove ultra-rare
    mask = (mean_abd > 0) & (occ > 0)
    mean_abd = mean_abd[mask]
    occ = occ[mask]

    m, rmse, pred = neutral_sloan_fit(mean_abd, occ, n_reads)

    # classify deviations (simple): observed - pred
    delta = occ - pred
    q = np.quantile(np.abs(delta), 0.95)
    cls = np.where(delta > q, "Above (selection/dispersion)", np.where(delta < -q, "Below (filtering)", "Neutral"))

    plot_df = pd.DataFrame(
        {
            "mean_abd": mean_abd,
            "occ": occ,
            "pred": pred,
            "class": cls,
        }
    )
    plot_df["log10_mean_abd"] = np.log10(plot_df["mean_abd"])

    fig, ax = plt.subplots(1, 1, figsize=(6, 4.5), constrained_layout=True)
    sns.scatterplot(
        data=plot_df.sample(min(20000, len(plot_df)), random_state=0),
        x="log10_mean_abd",
        y="occ",
        hue="class",
        palette={"Neutral": "#4C72B0", "Above (selection/dispersion)": "#C44E52", "Below (filtering)": "#8172B2"},
        s=10,
        linewidth=0,
        alpha=0.75,
        ax=ax,
    )
    # prediction curve
    xs = np.linspace(plot_df["log10_mean_abd"].min(), plot_df["log10_mean_abd"].max(), 200)
    meanx = 10**xs
    pred_curve = 1.0 - np.exp(-m * meanx * n_reads)
    ax.plot(xs, pred_curve, color="black", linewidth=1.5, label="Neutral fit")
    ax.set_title(f"Sloan-like neutral model (top features)\nfit m={m:.2e}, RMSE={rmse:.3f}")
    ax.set_xlabel("log10(mean relative abundance)")
    ax.set_ylabel("Occupancy (prevalence)")
    ax.legend(frameon=False, loc="lower right")

    for ext in ["pdf", "png"]:
        fig.savefig(outdir / f"Fig6_neutral_model_sloan_like_topFeatures.{ext}", bbox_inches="tight")
    plt.close(fig)


def fig_ddr_no_latlon(outdir: Path, pcs_tsv: Path, meta_tsv: Path, n_pairs: int = 200000):
    """
    Distance-decay surrogate without coordinates:
    compare within-country / between-country (and within-region / between-region) distance distributions in CLR-PC space.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    pcs = pd.read_csv(pcs_tsv, sep="\t")
    meta = pd.read_csv(meta_tsv, sep="\t")
    meta["sample_id"] = meta["sample_id"].astype(str)
    pcs["sample_id"] = pcs["sample_id"].astype(str)
    d = meta.merge(pcs, on="sample_id", how="inner")
    for c in d.columns:
        if c in {"sample_id"}:
            continue
        d[c] = d[c].map(norm_missing)

    pc_cols = [c for c in d.columns if c.startswith("PC")]
    X = d[pc_cols].to_numpy(dtype=np.float32)
    n = X.shape[0]
    rng = np.random.default_rng(0)

    def sample_same_diff_pairs(labels: np.ndarray, n_each: int):
        """Randomly sample pairs: n_each same-label pairs + n_each different-label pairs (both exclude NA)."""
        lab = labels.astype(object)
        same_i: list[int] = []
        same_j: list[int] = []
        diff_i: list[int] = []
        diff_j: list[int] = []
        tries = 0
        while (len(same_i) < n_each or len(diff_i) < n_each) and tries < max(2000, n_each * 50):
            batch = 8000
            i = rng.integers(0, n, size=batch)
            j = rng.integers(0, n, size=batch)
            ok = i != j
            i, j = i[ok], j[ok]
            li, lj = lab[i], lab[j]
            valid = (li != "NA") & (lj != "NA")
            i, j = i[valid], j[valid]
            li, lj = lab[i], lab[j]
            same_mask = li == lj
            diff_mask = li != lj
            if len(same_i) < n_each:
                si = i[same_mask]
                sj = j[same_mask]
                need = n_each - len(same_i)
                take = min(need, len(si))
                if take:
                    same_i.extend(si[:take].tolist())
                    same_j.extend(sj[:take].tolist())
            if len(diff_i) < n_each:
                di = i[diff_mask]
                dj = j[diff_mask]
                need = n_each - len(diff_i)
                take = min(need, len(di))
                if take:
                    diff_i.extend(di[:take].tolist())
                    diff_j.extend(dj[:take].tolist())
            tries += 1
        return np.asarray(same_i, dtype=int), np.asarray(same_j, dtype=int), np.asarray(diff_i, dtype=int), np.asarray(diff_j, dtype=int)

    def pair_dist(ii: np.ndarray, jj: np.ndarray) -> np.ndarray:
        return np.sqrt(((X[ii] - X[jj]) ** 2).sum(axis=1))

    n_each = max(2000, min(50000, n_pairs // 2))

    def plot_panel(ax, labels: np.ndarray, title: str, a: str, b: str):
        si, sj, di, dj = sample_same_diff_pairs(labels, n_each)
        ds = pair_dist(si, sj)
        dd = pair_dist(di, dj)
        plot = pd.DataFrame(
            {
                "distance": np.concatenate([ds, dd]),
                "group": [a] * len(ds) + [b] * len(dd),
            }
        )
        sns.violinplot(
            data=plot,
            x="group",
            y="distance",
            hue="group",
            inner="box",
            cut=0,
            ax=ax,
            palette=["#4C72B0", "#55A868"],
            legend=False,
        )
        ax.set_title(title)
        ax.set_xlabel("")
        ax.set_ylabel("Euclidean distance on PCs")
        for lab in ax.get_xticklabels():
            lab.set_rotation(15)
            lab.set_horizontalalignment("right")

    country = d["Country"].astype(str).to_numpy()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    plot_panel(axes[0], country, "Within-country vs Between-country", "Within-country", "Between-country")
    if "Region" in d.columns:
        region = d["Region"].astype(str).to_numpy()
        plot_panel(axes[1], region, "Within-Region vs Between-Region", "Within-Region", "Between-Region")
    else:
        axes[1].axis("off")
        axes[1].set_title("Region column missing")
    fig.suptitle("Distance-decay surrogate (CLR-PC space)", y=1.02)
    for ext in ["pdf", "png"]:
        fig.savefig(outdir / f"Fig7_ddr_within_between_country_region_violin.{ext}", bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, default=Path("figures_final"))
    ap.add_argument("--schemeA-combined", type=Path, default=Path("after_mmuphin_clr_adonis2_schemeA_combined.tsv"))
    ap.add_argument("--counts-top", type=Path, default=Path("pvca_Y_count.tsv"))
    ap.add_argument("--pcs", type=Path, default=Path("after_mmuphin_clr_pcs_k20.tsv"))
    ap.add_argument("--meta", type=Path, default=Path("pvca_sample_meta_full.tsv"))
    ap.add_argument("--ddr-pairs", type=int, default=200000)
    args = ap.parse_args()

    ensure_dir(args.outdir)
    set_pub_rcparams()

    fig_permanova_schemeA(args.outdir, args.schemeA_combined)
    fig_neutral_model_top30000(args.outdir, args.counts_top)
    fig_ddr_no_latlon(args.outdir, args.pcs, args.meta, n_pairs=args.ddr_pairs)

    print(f"[DONE] wrote part2 figures to {args.outdir.resolve()}")


if __name__ == "__main__":
    main()

