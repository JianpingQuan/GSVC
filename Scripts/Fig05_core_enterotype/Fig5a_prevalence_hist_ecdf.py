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
"""
Core / satellite: prevalence distribution + core≥0.5 / ≥0.9 vOTU counts;
Virus annotation (phylum/class/order) and host prediction (phylum/family/genus) pie charts (PDF+PNG, 300 dpi).

Data:
  - FINAL_all_projects_Count_matrix.tsv
  - high_medium_quality_viral_rep_seq_virus_summary.tsv (taxonomy; seq_name = matrix Contig)
  - high/medium_quality_Host_prediction_to_genome_m90.csv (best host per vOTU by Confidence)

Output: 06_core_virome_figs/figures/ and tables/
Environment variable PHAGE_PROJECT_ROOT
"""
from __future__ import annotations

import csv
import os
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

_THIS = Path(__file__).resolve()
_DEFAULT_ROOT = _THIS.parent.parent
ROOT = Path(os.environ.get("PHAGE_PROJECT_ROOT", str(_DEFAULT_ROOT))).expanduser().resolve()
if not ROOT.exists():
    ROOT = _DEFAULT_ROOT

OUT = _THIS.parent
FIG = OUT / "figures"
TAB = OUT / "tables"
for d in (FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

MATRIX = ROOT / "FINAL_all_projects_Count_matrix.tsv"
VIRUS_SUMMARY = ROOT / "high_medium_quality_viral_rep_seq_virus_summary.tsv"
HOST_HQ = ROOT / "high_quality_Host_prediction_to_genome_m90.csv"
HOST_MQ = ROOT / "medium_quality_Host_prediction_to_genome_m90.csv"

CHUNK_ROWS = 2500
DETECT_GT = 0.0
TOP_PIE = 10
PREV_THRESHOLDS = (0.5, 0.9)

# Nature / general journal friendly: low saturation, distinguishable
NATURE_COLORS = [
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
    "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
    "#499894", "#E377C2", "#8C564B", "#BCBD22", "#17BECF",
    "#9467BD", "#7F7F7F", "#D62728", "#1F77B4", "#2CA02C",
]


def apply_pub_style() -> None:
    import matplotlib as mpl

    mpl.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.dpi": 300,
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "Liberation Sans"],
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "axes.titleweight": "normal",
            "axes.linewidth": 0.8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.color": "#E8E8E8",
            "grid.linewidth": 0.5,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def parse_sample_ids(header: list[str]) -> list[str]:
    ids = []
    for c in header[1:]:
        c = str(c).strip()
        if c.endswith("_Count"):
            ids.append(c[: -len("_Count")])
        else:
            ids.append(c)
    return ids


def stream_prevalence(matrix_path: Path) -> tuple[list[str], np.ndarray, int]:
    """Return (contig_ids, prevalence vector, n_samples)."""
    with matrix_path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        header = f.readline().rstrip("\n\r").split("\t")
    sample_ids = parse_sample_ids(header)
    n_s = len(sample_ids)
    if n_s == 0:
        raise ValueError("No sample column")

    contigs: list[str] = []
    prevs: list[float] = []

    with matrix_path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader, None)
        buf: list[list[str]] = []
        for row in reader:
            if not row:
                continue
            buf.append(row)
            if len(buf) >= CHUNK_ROWS:
                _flush_chunk(buf, sample_ids, n_s, contigs, prevs)
                buf = []
        if buf:
            _flush_chunk(buf, sample_ids, n_s, contigs, prevs)

    return contigs, np.array(prevs, dtype=np.float64), n_s


def _flush_chunk(
    buf: list[list[str]],
    sample_ids: list[str],
    n_s: int,
    contigs: list[str],
    prevs: list[float],
) -> None:
    for row in buf:
        votu = row[0].strip()
        vals = np.fromstring("\t".join(row[1:]), sep="\t", dtype=np.float64)
        if vals.size < n_s:
            vals = np.pad(vals, (0, n_s - vals.size), constant_values=0.0)
        elif vals.size > n_s:
            vals = vals[:n_s]
        prev = float(np.mean(vals > DETECT_GT))
        contigs.append(votu)
        prevs.append(prev)


def parse_virus_ranks(tax: str) -> dict[str, str]:
    out = {"phylum": "Unclassified", "class": "Unclassified", "order": "Unclassified"}
    if pd.isna(tax) or not str(tax).strip():
        return out
    parts = [x.strip() for x in str(tax).split(";") if x.strip()]
    for p in parts:
        if re.search(r"viricota$", p, re.I):
            out["phylum"] = p
        if re.search(r"viricetes$", p, re.I):
            out["class"] = p
        if re.search(r"virales$", p, re.I):
            out["order"] = p
    return out


def parse_host_rank(tax: str, rank: str) -> str:
    if pd.isna(tax) or not str(tax).strip():
        return "Unclassified"
    pref = {"phylum": "p__", "family": "f__", "genus": "g__"}[rank]
    for part in str(tax).split(";"):
        p = part.strip()
        if p.startswith(pref):
            s = p[len(pref) :].strip()
            return s if s else "Unclassified"
    return "Unclassified"


def load_virus_taxonomy_map() -> dict[str, dict[str, str]]:
    """seq_name -> {phylum, class, order}"""
    out: dict[str, dict[str, str]] = {}
    for ch in pd.read_csv(VIRUS_SUMMARY, sep="\t", chunksize=80_000, usecols=["seq_name", "taxonomy"]):
        for sid, tax in zip(
            ch["seq_name"].astype(str).str.strip(),
            ch["taxonomy"],
        ):
            out[sid] = parse_virus_ranks(tax)
    return out


def build_virus_best_host() -> dict[str, str]:
    """Virus contig -> Host taxonomy string (max confidence)."""
    best: dict[str, tuple[float, str]] = {}
    cols = ["Virus", "Host taxonomy", "Confidence score"]
    for path in (HOST_HQ, HOST_MQ):
        p = Path(path)
        if not p.exists():
            continue
        for ch in pd.read_csv(p, chunksize=400_000, usecols=lambda c: c in cols, low_memory=False):
            ch["Confidence score"] = pd.to_numeric(ch["Confidence score"], errors="coerce").fillna(0.0)
            for i in range(len(ch)):
                v = str(ch["Virus"].iloc[i]).strip()
                if not v or v == "nan":
                    continue
                conf = float(ch["Confidence score"].iloc[i])
                tax = ch["Host taxonomy"].iloc[i]
                old = best.get(v)
                if old is None or conf > old[0]:
                    best[v] = (conf, str(tax) if pd.notna(tax) else "")
    return {v: t[1] for v, t in best.items()}


def counter_for_votus(
    votu_ids: list[str],
    virus_map: dict[str, dict[str, str]],
    host_map: dict[str, str],
    virus_key: str,
    host_rank: str | None,
) -> Counter:
    c: Counter = Counter()
    for vid in votu_ids:
        if host_rank is None:
            rk = virus_map.get(vid, {}).get(virus_key, "Unclassified")
            if not rk or rk == "":
                rk = "Unclassified"
            c[rk] += 1
        else:
            tax = host_map.get(vid, "")
            rk = parse_host_rank(tax, host_rank)
            c[rk] += 1
    return c


def collapse_counter(ct: Counter, top_n: int) -> tuple[list[str], np.ndarray]:
    if not ct:
        return ["Unclassified"], np.array([1.0])
    items = ct.most_common()
    if len(items) <= top_n:
        labs = [x[0] for x in items]
        vals = np.array([x[1] for x in items], dtype=float)
    else:
        top = items[:top_n]
        other = sum(x[1] for x in items[top_n:])
        labs = [x[0] for x in top] + ["Other"]
        vals = np.array([x[1] for x in top] + [other], dtype=float)
    vals = vals / vals.sum()
    # shorten labels
    labs = [str(l)[:42] + "…" if len(str(l)) > 42 else str(l) for l in labs]
    return labs, vals


def plot_pie_ax(ax, labels: list[str], fracs: np.ndarray, title: str) -> None:
    n = len(labels)
    colors = [NATURE_COLORS[i % len(NATURE_COLORS)] for i in range(n)]
    wedges, _texts, autotexts = ax.pie(
        fracs,
        labels=None,
        autopct=lambda pct: f"{pct:.1f}%" if pct >= 3 else "",
        colors=colors,
        pctdistance=0.72,
        wedgeprops={"linewidth": 0.6, "edgecolor": "white"},
        textprops={"size": 8, "color": "#222222"},
    )
    for t in autotexts:
        t.set_fontsize(7)
    ax.set_title(title, fontsize=10, pad=8)
    ax.legend(
        wedges,
        [f"{l} ({f*100:.1f}%)" for l, f in zip(labels, fracs)],
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=7,
        frameon=False,
        borderaxespad=0.0,
    )


def main() -> None:
    import argparse
    import matplotlib

    ap = argparse.ArgumentParser(description="Core virome prevalence + taxonomy/host pies")
    ap.add_argument(
        "--skip-full-prev-table",
        action="store_true"
    )
    args = ap.parse_args()

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    apply_pub_style()

    if not MATRIX.exists():
        raise FileNotFoundError(MATRIX)
    if not VIRUS_SUMMARY.exists():
        raise FileNotFoundError(VIRUS_SUMMARY)

    print("prevalence …", flush=True)
    contigs, prev, n_samples = stream_prevalence(MATRIX)
    print("vOTU", len(contigs), "Sample:", n_samples, flush=True)

    if not args.skip_full_prev_table:
        pd.DataFrame({"vOTU": contigs, "prevalence": prev}).to_csv(
            TAB / "all_vOTU_prevalence_full.tsv", sep="\t", index=False
        )

    n50 = int(np.sum(prev >= 0.5))
    n90 = int(np.sum(prev >= 0.9))
    with open(TAB / "core_prevalence_counts_summary.txt", "w", encoding="utf-8") as f:
        f.write(f"n_samples\t{n_samples}\n")
        f.write(f"n_vOTU_total\t{len(contigs)}\n")
        f.write(f"n_vOTU_prevalence_ge_0.5\t{n50}\n")
        f.write(f"n_vOTU_prevalence_ge_0.9\t{n90}\n")

    set50 = [contigs[i] for i in range(len(contigs)) if prev[i] >= 0.5]
    set90 = [contigs[i] for i in range(len(contigs)) if prev[i] >= 0.9]
    pd.DataFrame({"vOTU": set50}).to_csv(TAB / "core_vOTU_prev_ge0.5.tsv", sep="\t", index=False)
    pd.DataFrame({"vOTU": set90}).to_csv(TAB / "core_vOTU_prev_ge0.9.tsv", sep="\t", index=False)

    # Fig: prevalence histogram + ECDF (publication line weights and margins)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), constrained_layout=True)
    axes[0].hist(prev, bins=50, range=(0, 1), color="#4E79A7", edgecolor="white", linewidth=0.35)
    axes[0].axvline(0.5, color="#C44E52", linestyle="--", linewidth=1.2, label="prevalence = 0.5")
    axes[0].axvline(0.9, color="#55A868", linestyle="--", linewidth=1.2, label="prevalence = 0.9")
    axes[0].set_xlabel("Prevalence (fraction of samples)")
    axes[0].set_ylabel("Number of vOTUs")
    axes[0].set_title("Prevalence distribution (all vOTUs)")
    axes[0].legend(frameon=False, loc="upper right")
    axes[0].set_axisbelow(True)
    pv = np.sort(prev)
    y = np.arange(1, len(pv) + 1) / len(pv)
    axes[1].plot(pv, y, color="#4E79A7", linewidth=1.4)
    axes[1].axvline(0.5, color="#C44E52", linestyle="--", linewidth=1.2)
    axes[1].axvline(0.9, color="#55A868", linestyle="--", linewidth=1.2)
    axes[1].set_xlabel("Prevalence")
    axes[1].set_ylabel("ECDF")
    axes[1].set_title("Empirical cumulative distribution")
    axes[1].set_xlim(0, 1)
    axes[1].set_ylim(0, 1.02)
    axes[1].set_axisbelow(True)
    fig.suptitle(
        f"vOTU prevalence across {n_samples} samples (detected if count > {DETECT_GT})",
        fontsize=11,
        y=1.03,
    )
    for ext in (".pdf", ".png"):
        fig.savefig(
            FIG / f"Fig_prevalence_hist_ecdf{ext}",
            dpi=300 if ext == ".png" else None,
            bbox_inches="tight",
            pad_inches=0.03,
        )
    plt.close(fig)

    print(" taxonomy …", flush=True)
    virus_map = load_virus_taxonomy_map()
    print("host prediction …", flush=True)
    host_map = build_virus_best_host()

    def make_pie_panel(
        suffix: str,
        votu_list: list[str],
        label: str,
    ) -> None:
        fig, axes = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=False)
        # row0 virus: phylum, class, order
        for j, (vk, tit) in enumerate(
            [("phylum", "Virus phylum"), ("class", "Virus class"), ("order", "Virus order")]
        ):
            ct = counter_for_votus(votu_list, virus_map, host_map, vk, None)
            labs, fr = collapse_counter(ct, TOP_PIE)
            plot_pie_ax(axes[0, j], labs, fr, f"{label}\n{tit}\n(n={len(votu_list)} vOTUs)")
        # row1 host: phylum, family, genus
        for j, (hr, tit) in enumerate(
            [("phylum", "Host phylum"), ("family", "Host family"), ("genus", "Host genus")]
        ):
            ct = counter_for_votus(votu_list, virus_map, host_map, "", hr)
            labs, fr = collapse_counter(ct, TOP_PIE)
            plot_pie_ax(axes[1, j], labs, fr, f"{label}\n{tit}")

        fig.suptitle(
            f"{label}: composition by vOTU count (top {TOP_PIE} categories + Other)",
            fontsize=12,
            y=0.995,
        )
        fig.subplots_adjust(left=0.06, right=0.82, top=0.90, bottom=0.06, wspace=0.35, hspace=0.22)
        for ext in (".pdf", ".png"):
            fig.savefig(
                FIG / f"Fig_core_{suffix}_virus_host_taxonomy_pies{ext}",
                dpi=300 if ext == ".png" else None,
                bbox_inches="tight",
                pad_inches=0.05,
            )
        plt.close(fig)

    make_pie_panel("prev_ge0.5", set50, "Core-like (prevalence ≥ 0.5)")
    make_pie_panel("prev_ge0.9", set90, "Core (prevalence ≥ 0.9)")

    # long table composition
    rows = []
    for tag, votu_list in [("ge0.5", set50), ("ge0.9", set90)]:
        for vk in ("phylum", "class", "order"):
            ct = counter_for_votus(votu_list, virus_map, host_map, vk, None)
            tot = sum(ct.values()) or 1
            for name, c in ct.most_common():
                rows.append(
                    {"threshold": tag, "type": "virus", "rank": vk, "label": name, "n_vOTU": c, "fraction": c / tot}
                )
        for hr in ("phylum", "family", "genus"):
            ct = counter_for_votus(votu_list, virus_map, host_map, "", hr)
            tot = sum(ct.values()) or 1
            for name, c in ct.most_common():
                rows.append(
                    {"threshold": tag, "type": "host", "rank": hr, "label": name, "n_vOTU": c, "fraction": c / tot}
                )
    pd.DataFrame(rows).to_csv(TAB / "core_taxonomy_host_composition_long.tsv", sep="\t", index=False)

    print("Finished:", FIG, TAB, flush=True)


if __name__ == "__main__":
    main()
