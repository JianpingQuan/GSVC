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
"""Summarize PRJNA1010706 cohort meta tables as publication-style figures."""
from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np

BASE = Path(__file__).resolve().parent
FIG_DIR = BASE / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

MERGED = BASE / "PRJNA1010706_pheno_meta_merged.tsv"
SUM_AGE = BASE / "summary_n_by_cohort_age.tsv"

# Matches manuscript arms (GG + G1–G6) ↔ Treatment_x in meta / PigID prefix
ARM_ORDER = ["GG", "G1", "G2", "G3", "G4", "G5", "G6"]
TREATMENT_TO_ARM = {
    "Trimethoprim/Sulfamethoxazole": "G1",
    "Colistin": "G2",
    "Oral attenuated vaccine": "G3",
    "Gentamicin": "G4",
    "Control with water acidification": "G5",
    "Control with untreated water": "G6",
    "Amoxicillin": "GG",
}
ARM_TO_TREATMENT = {arm: tx for tx, arm in TREATMENT_TO_ARM.items()}



def treatments_in_arm_order() -> list[str]:
    return [ARM_TO_TREATMENT[a] for a in ARM_ORDER]


def parse_pheno_float(raw: str | None) -> float | None:
    x = (raw or "").strip()
    if x in ("", "/", "NA", "na"):
        return None
    try:
        return float(x)
    except ValueError:
        return None


def dedupe_one_row_per_pig(merged_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Phenotypes are identical across SRRs from the same PigID; stats deduplicate to one row per pig."""
    seen: set[str] = set()
    out: list[dict[str, str]] = []
    for r in merged_rows:
        pid = (r.get("cohort_PigID") or r.get("PigID") or "").strip()
        if not pid or pid in seen:
            continue
        seen.add(pid)
        out.append(r)
    return out


def _set_cn_font() -> None:
    candidates = ["Microsoft YaHei", "SimHei", "PingFang SC", "Noto Sans CJK SC", "Arial Unicode MS"]
    available = {f.name for f in fm.fontManager.ttflist}
    for name in candidates:
        if name in available:
            plt.rcParams["font.sans-serif"] = [name, "DejaVu Sans"]
            break
    else:
        plt.rcParams["font.sans-serif"] = ["DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def save(fig: plt.Figure, stem: str) -> None:
    for ext, kw in ((".png", {"dpi": 320}), (".pdf", {})):
        p = FIG_DIR / f"{stem}{ext}"
        fig.savefig(p, bbox_inches="tight", **kw)
    plt.close(fig)


def plot_sample_counts_by_arm(merged_rows: list[dict[str, str]], stem: str = "Fig1_sample_counts_by_treatment") -> None:
    c = Counter((r.get("paper_arm") or "").strip() for r in merged_rows if (r.get("paper_arm") or "").strip())
    counts = [int(c.get(a, 0)) for a in ARM_ORDER]
    labels = [ARM_LABEL_ZH[a] for a in ARM_ORDER]
    fig, ax = plt.subplots(figsize=(9.0, 5.8))
    y = np.arange(len(ARM_ORDER))
    ax.barh(y, counts, color="#2E6F95", edgecolor="white", linewidth=0.6)
    ax.set_yticks(y, labels, fontsize=10)
    ax.set_xlabel("SRR", fontsize=12)
    ax.set_title(
        "PRJNA1010706",
        fontsize=13,
        fontweight="600",
    )
    for yi, n in zip(y, counts):
        ax.text(n + 0.8, yi, str(n), va="center", fontsize=10, color="#333")
    ax.set_xlim(0, max(counts) * 1.18 if counts else 1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    save(fig, stem)


def plot_age(rows: list[dict[str, str]], stem: str = "Fig2_sample_counts_by_age") -> None:
    ages = [r["cohort_Age"] for r in rows if r.get("cohort_Age") and r["cohort_Age"] != "(empty)"]
    counts = [int(r["n_rows_in_pheno_file"]) for r in rows if r.get("cohort_Age") and r["cohort_Age"] != "(empty)"]
    pairs = sorted(zip(ages, counts), key=lambda x: int(x[0]) if x[0].isdigit() else 999)
    ages, counts = [p[0] for p in pairs], [p[1] for p in pairs]
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    x = np.arange(len(ages))
    bars = ax.bar(x, counts, color="#6B4E71", edgecolor="white", linewidth=0.6)
    ax.set_xticks(x, [f"{a} d" for a in ages], fontsize=11)
    ax.set_ylabel("SRR", fontsize=12)
    ax.set_title("PRJNA1010706：Age distribution（cohort table）", fontsize=13, fontweight="600")
    for b, c in zip(bars, counts):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 1.2, str(c), ha="center", fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.subplots_adjust(bottom=0.28)
    fig.text(
        0.5,
        0.02,
        fontsize=8.2,
        ha="center",
        color="#333",
    )
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    save(fig, stem)


def plot_treatment_age_heatmap(merged_rows: list[dict[str, str]], stem: str = "Fig3_treatment_by_age_counts") -> None:
    ages = sorted(
        {r.get("cohort_Age", "").strip() for r in merged_rows if r.get("cohort_Age", "").strip()},
        key=lambda a: int(a) if a.isdigit() else 999,
    )
    mat = np.zeros((len(ARM_ORDER), len(ages)), dtype=int)
    ai = {a: j for j, a in enumerate(ages)}
    gi = {a: i for i, a in enumerate(ARM_ORDER)}
    for r in merged_rows:
        arm = (r.get("paper_arm") or "").strip()
        a = (r.get("cohort_Age") or "").strip()
        if arm in gi and a in ai:
            mat[gi[arm], ai[a]] += 1
    fig, ax = plt.subplots(figsize=(7.4, 6.4))
    im = ax.imshow(mat, aspect="auto", cmap="YlOrBr", vmin=0, vmax=mat.max() or 1)
    ax.set_xticks(np.arange(len(ages)), [f"{a} d" for a in ages], fontsize=10)
    ax.set_yticks(np.arange(len(ARM_ORDER)), [ARM_LABEL_ZH[a] for a in ARM_ORDER], fontsize=9)
    ax.set_xlabel("Age（cohort table）", fontsize=12)
    ax.set_ylabel("Arm（Fig.1）", fontsize=12)
    ax.set_title("PRJNA1010706：Arm × Age", fontsize=13, fontweight="600")
    vmax = float(mat.max()) if mat.size else 1.0
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat[i, j]
            color = "white" if v > vmax / 2 else "#222"
            ax.text(j, i, str(int(v)), ha="center", va="center", fontsize=10, color=color)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Sample number", fontsize=11)
    fig.tight_layout()
    save(fig, stem)


def plot_phenotype_by_treatment(
    pig_rows: list[dict[str, str]], stem: str = "Fig4_phenotype_by_treatment"
) -> None:
    treat_order = treatments_in_arm_order()
    xlabels = [ARM_LABEL_ZH[TREATMENT_TO_ARM[tr]] for tr in treat_order]
    metrics = [
        ("cohort_Initial weight", "#1B6B5E"),
        ("cohort_Final weight", "#3D8B72"),
        ("cohort_ADWG", "ADWG (kg/d)", "#6BB39B"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 5.8))
    for ax, (col, ylabel, color) in zip(axes, metrics):
        data, pos = [], []
        for j, tr in enumerate(treat_order):
            vals = [
                parse_pheno_float(r.get(col))
                for r in pig_rows
                if (r.get("Treatment_x") or "").strip() == tr
            ]
            vals = [v for v in vals if v is not None]
            if vals:
                data.append(vals)
                pos.append(j)
        if not data:
            ax.set_visible(False)
            continue
        bp = ax.boxplot(
            data,
            positions=pos,
            widths=0.55,
            patch_artist=True,
            medianprops={"color": "#222", "linewidth": 1.4},
            boxprops={"facecolor": color, "edgecolor": "#333", "linewidth": 0.8},
            whiskerprops={"color": "#333"},
            capprops={"color": "#333"},
            flierprops={"marker": "o", "markersize": 4, "alpha": 0.55},
        )
        for b in bp["boxes"]:
            b.set(alpha=0.82)
        ax.set_xticks(range(len(treat_order)))
        ax.set_xticklabels(xlabels, rotation=28, ha="right", fontsize=8.5)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", linestyle=":", alpha=0.45)
    )
    fig.tight_layout()
    save(fig, stem)


def plot_phenotype_by_age_srr(
    merged_rows: list[dict[str, str]], stem: str = "Fig5_phenotype_by_age_SRR"
) -> None:
    """Plot age and phenotypes at SRR level; phenotype values repeat when the same pig has multiple time points."""
    ages = sorted(
        {r.get("cohort_Age", "").strip() for r in merged_rows if r.get("cohort_Age", "").strip()},
        key=lambda a: int(a) if a.isdigit() else 999,
    )
    metrics = [
        ("cohort_Initial weight",  "#4A5899"),
        ("cohort_Final weight", "#6B7BC4"),
        ("cohort_ADWG", "ADWG (kg/d)", "#8FA0D9"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 5.2))
    x = np.arange(len(ages))
    for ax, (col, ylabel, color) in zip(axes, metrics):
        parts = []
        for a in ages:
            vals = [
                parse_pheno_float(r.get(col))
                for r in merged_rows
                if (r.get("cohort_Age") or "").strip() == a
            ]
            parts.append([v for v in vals if v is not None])
        bp = ax.boxplot(
            parts,
            positions=x,
            widths=0.5,
            patch_artist=True,
            medianprops={"color": "#222", "linewidth": 1.4},
            boxprops={"facecolor": color, "edgecolor": "#333", "linewidth": 0.8},
            whiskerprops={"color": "#333"},
            capprops={"color": "#333"},
            flierprops={"marker": "o", "markersize": 3, "alpha": 0.45},
        )
        for b in bp["boxes"]:
            b.set(alpha=0.82)
        ax.set_xticks(x, [f"{a} d" for a in ages], fontsize=11)
        ax.set_xlabel("Age（cohort table）", fontsize=11)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", linestyle=":", alpha=0.45)
    fig.tight_layout()
    save(fig, stem)


def plot_initial_vs_final_scatter(
    pig_rows: list[dict[str, str]], stem: str = "Fig6_initial_vs_final_by_treatment"
) -> None:
    treat_order = treatments_in_arm_order()
    cmap = plt.matplotlib.colormaps.get_cmap("tab10")
    colors = cmap(np.linspace(0, 0.92, max(len(treat_order), 1), endpoint=True))
    fig, ax = plt.subplots(figsize=(7.6, 6.4))
    for i, tr in enumerate(treat_order):
        arm = TREATMENT_TO_ARM[tr]
        xs, ys = [], []
        for r in pig_rows:
            if (r.get("Treatment_x") or "").strip() != tr:
                continue
            ix = parse_pheno_float(r.get("cohort_Initial weight"))
            iy = parse_pheno_float(r.get("cohort_Final weight"))
            if ix is not None and iy is not None:
                xs.append(ix)
                ys.append(iy)
        if not xs:
            continue
        ax.scatter(
            xs,
            ys,
            s=52,
            alpha=0.78,
            color=colors[i],
            edgecolors="white",
            linewidths=0.6,
            label=ARM_LABEL_ZH[arm].replace("\n", " "),
        )
    ax.legend(loc="upper left", fontsize=7.0, frameon=True, framealpha=0.92)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(linestyle=":", alpha=0.45)
    fig.tight_layout()
    save(fig, stem)


def main() -> None:
    _set_cn_font()
    if not MERGED.exists():
        raise SystemExit(f"Missing merged table: {MERGED}")

    merged = read_tsv(MERGED)
    if not any((r.get("paper_arm") or "").strip() for r in merged):
        raise SystemExit("Merged table has no paper_arm column; run build_merged_with_paper_arm.py first.")

    plot_sample_counts_by_arm(merged)
    plot_age(read_tsv(SUM_AGE))
    plot_treatment_age_heatmap(merged)

    pigs = dedupe_one_row_per_pig(merged)
    plot_phenotype_by_treatment(pigs)
    plot_phenotype_by_age_srr(merged)
    plot_initial_vs_final_scatter(pigs)

    print("Wrote figures to", FIG_DIR)


if __name__ == "__main__":
    main()
