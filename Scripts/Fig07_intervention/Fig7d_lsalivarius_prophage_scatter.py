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


# -*- coding: utf-8 -*-
"""Seven scatter plots: linear fit + gray 95% CI; Science/AAAS-style; Spearman rho & p only."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import spearmanr
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import config  # noqa: E402

FIG = config.FIGURES / "lsalivarius_prophage_evidence"
MASTER = config.RESULTS / "lsalivarius_prophage_evidence" / "pig_master_evidence.tsv"

LSAL_COL = "log1p_Ligilactobacillus_salivarius"
MIN_Lsal = 0.2
MIN_PREV = 0.2
SUBSET_NAME = "Lsal_detected_log1p_gt0.2"

# AAAS / Science-inspired (colorblind-friendly Wong-like + AAAS accents)
SCIENCE = {
    "point": "#0072B2",
    "point_edge": "#FFFFFF",
    "line": "#D55E00",
    "ci": "#B8B8B8",
    "text": "#1A1A1A",
    "spine": "#333333",
}

PLOT_SPECS = [
    (LSAL_COL, "ADWG", "scatter_Lsalivarius_vs_ADWG_Lsal_detected_log1p_gt0", False),
    ("CRR456631_viral_contig_2066", "ADWG", "scatter_CRR456631_viral_contig_2066_vs_ADWG_Lsal_detected_log1p_gt0", True),
    ("CRR456631_viral_contig_2066", LSAL_COL, "scatter_CRR456631_viral_contig_2066_vs_Lsal_Lsal_detected_log1p_gt0", True),
    ("SRR15732359_viral_contig_2826", "ADWG", "scatter_SRR15732359_viral_contig_2826_vs_ADWG_Lsal_detected_log1p_gt0", True),
    ("SRR15732359_viral_contig_2826", LSAL_COL, "scatter_SRR15732359_viral_contig_2826_vs_Lsal_Lsal_detected_log1p_gt0", True),
    ("SRR31546707_viral_contig_21744", "ADWG", "scatter_SRR31546707_viral_contig_21744_vs_ADWG_Lsal_detected_log1p_gt0", True),
    ("SRR31546707_viral_contig_21744", LSAL_COL, "scatter_SRR31546707_viral_contig_21744_vs_Lsal_Lsal_detected_log1p_gt0", True),
]


def apply_science_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 9,
            "axes.labelsize": 10,
            "axes.titlesize": 10,
            "axes.titleweight": "normal",
            "axes.labelcolor": SCIENCE["text"],
            "axes.edgecolor": SCIENCE["spine"],
            "axes.linewidth": 0.8,
            "xtick.color": SCIENCE["spine"],
            "ytick.color": SCIENCE["spine"],
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def spearman_rho_p(sub: pd.DataFrame, x_col: str, y_col: str) -> tuple[float, float]:
    rho, p = spearmanr(sub[x_col], sub[y_col])
    return float(rho), float(p)


def format_rho_p(rho: float, p: float) -> str:
    if not np.isfinite(rho) or not np.isfinite(p):
        return "ρ = NA, p = NA"
    ps = "p < 0.001" if p < 0.001 else f"p = {p:.4f}"
    return f"ρ = {rho:.3f}, {ps}"


def add_linear_ci(ax, x: np.ndarray, y: np.ndarray) -> None:
    if len(x) < 3:
        return
    fit = sm.OLS(y, sm.add_constant(x)).fit()
    x_grid = np.linspace(float(np.min(x)), float(np.max(x)), 100)
    sf = fit.get_prediction(sm.add_constant(x_grid)).summary_frame(alpha=0.05)
    ax.fill_between(
        x_grid,
        sf["mean_ci_lower"],
        sf["mean_ci_upper"],
        color=SCIENCE["ci"],
        alpha=0.45,
        linewidth=0,
        zorder=2,
    )
    ax.plot(x_grid, sf["mean"], color=SCIENCE["line"], lw=1.8, alpha=0.95, zorder=4)


def labels_for(x_col: str, y_col: str) -> tuple[str, str, str]:
    if x_col == LSAL_COL:
        xlab = "L. salivarius mean log1p TPM"
    else:
        xlab = f"{x_col} mean log1p TPM"
    ylab = "ADWG" if y_col == "ADWG" else "L. salivarius mean log1p TPM"
    if x_col == LSAL_COL and y_col == "ADWG":
        title = f"L. salivarius vs ADWG\n({SUBSET_NAME})"
    elif y_col == "ADWG":
        title = f"{x_col}\nvs ADWG ({SUBSET_NAME})"
    else:
        title = f"{x_col}\nvs L. salivarius ({SUBSET_NAME})"
    return xlab, ylab, title


def scatter_plot(sub: pd.DataFrame, x_col: str, y_col: str, out_stem: str) -> None:
    apply_science_style()
    sub = sub[[x_col, y_col]].dropna()
    sub = sub[np.isfinite(sub[x_col]) & np.isfinite(sub[y_col])]
    rho, p = spearman_rho_p(sub, x_col, y_col)
    xlab, ylab, title = labels_for(x_col, y_col)

    fig, ax = plt.subplots(figsize=(4.2, 4.0), dpi=150)
    x = sub[x_col].to_numpy()
    y = sub[y_col].to_numpy()

    ax.scatter(
        x,
        y,
        s=40,
        alpha=0.88,
        edgecolors=SCIENCE["point_edge"],
        linewidths=0.5,
        c=SCIENCE["point"],
        zorder=3,
    )
    add_linear_ci(ax, x, y)

    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title, fontsize=10, color=SCIENCE["text"])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.text(
        0.03,
        0.97,
        format_rho_p(rho, p),
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        color=SCIENCE["text"],
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white", alpha=0.92, edgecolor="#DDDDDD", lw=0.6),
    )

    fig.tight_layout()
    out = FIG / out_stem
    fig.savefig(out.with_suffix(".png"), dpi=300)
    try:
        fig.savefig(out.with_suffix(".pdf"))
    except OSError as e:
        print(f"  (PDF skipped — close file if open: {e})")
    plt.close(fig)
    print(f"Wrote {out_stem}.png/pdf  {format_rho_p(rho, p)}")


def main() -> None:
    if not MASTER.is_file():
        raise SystemExit(f"Run 25_lsalivarius_prophage_evidence_no_mediation.py first; missing {MASTER}")

    FIG.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(MASTER, sep="\t")
    base = df[df[LSAL_COL] > MIN_Lsal].copy()

    for x_col, y_col, stem, filt_x in PLOT_SPECS:
        sub = base
        if filt_x and x_col != LSAL_COL:
            sub = sub[sub[x_col] > MIN_PREV].copy()
        scatter_plot(sub, x_col, y_col, stem)

    print("Done — 7 scatter plots in", FIG)


if __name__ == "__main__":
    main()
