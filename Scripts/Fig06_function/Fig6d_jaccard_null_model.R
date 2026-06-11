# GSVC_github publication script
# =============================================================================
# GSVC paper figure script (standalone copy for publication reproducibility)
# Regenerate: python paper_figure_scripts/generate_paper_scripts.py
# =============================================================================
args <- commandArgs(trailingOnly = FALSE)
sp <- sub("^--file=", "", args[grep("^--file=", args)])
PHAGE_ROOT <- if (length(sp) && nzchar(sp)) {
  normalizePath(file.path(dirname(sp), "../.."), winslash = "/", mustWork = FALSE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
if (!file.exists(file.path(PHAGE_ROOT, "FINAL_all_projects_Count_matrix.tsv"))) {
  ev <- Sys.getenv("PHAGE_PROJECT_ROOT", "")
  if (nzchar(ev)) PHAGE_ROOT <- normalizePath(ev, winslash = "/")
}
OUT_DIR <- file.path(PHAGE_ROOT, "Scripts", "output")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(PHAGE_ROOT)


from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.patches import Circle, PathPatch, Wedge
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
import networkx as nx
import numpy as np
import pandas as pd


ROOT = (Path(os.environ.get('GSVC_ROOT') or Path(os.environ.get('PHAGE_PROJECT_ROOT') or Path(__file__).resolve().parents[2])).resolve() / 'Functional')
IN_DIR = ROOT / "summary_tables" / "host_phage_ko_coupling"
OUT_DIR = ROOT / "figures" / "host_phage_ko_coupling"

TAX_PATH = ROOT / "taxanomy_final.txt"
VIB_INDIV = ROOT / "VIBRANT_AMG_individuals_high_quality_viral_rep_seq.tsv"
CRISPR_BEST = ROOT / "summary_tables" / "amg_spacer_host" / "spacer_best_hit_per_spacer.tsv"

PALETTE_HEX = ["#8FB4DC", "#FFDD8E", "#70CDBE", "#AC99D2", "#7AC3DF", "#F5AA61", "#EB756C"]


def _set_plot_style() -> None:
    mpl.rcParams.update(
        {
            "figure.dpi": 140,
            "savefig.dpi": 300,
            "font.size": 11,
            "font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
            "axes.unicode_minus": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 1.0,
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "xtick.direction": "out",
            "ytick.direction": "out",
        }
    )


def safe_log2(x: float) -> float:
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return float("nan")
    if x <= 0:
        return float("nan")
    return math.log2(x)


@dataclass(frozen=True)
class OrCI:
    log2_or: float
    log2_lo: float
    log2_hi: float


def or_ci_haldane(a: int, b: int, c: int, d: int) -> OrCI:
    """
    2x2 odds ratio CI via log(OR) +/- 1.96*SE with Haldane-Anscombe correction.

    Table:
        phage+  phage-
    host+   a      b
    host-   c      d
    """
    aa = a + 0.5
    bb = b + 0.5
    cc = c + 0.5
    dd = d + 0.5
    or_hat = (aa * dd) / (bb * cc)
    se = math.sqrt(1.0 / aa + 1.0 / bb + 1.0 / cc + 1.0 / dd)
    lo = math.exp(math.log(or_hat) - 1.96 * se)
    hi = math.exp(math.log(or_hat) + 1.96 * se)
    return OrCI(log2_or=safe_log2(or_hat), log2_lo=safe_log2(lo), log2_hi=safe_log2(hi))


def fig_x_forest_directionality(path: Path, out_png: Path, out_pdf: Path) -> None:
    df = pd.read_csv(path, sep="\t", dtype={"KO": str})
    if df.empty:
        return

    # keep rows with numeric FDR
    df["q_bh_host_ge"] = pd.to_numeric(df["q_bh_host_ge"], errors="coerce")

    # compute CI on log2 scale, robust to inf OR
    cis: list[OrCI] = []
    for _, r in df.iterrows():
        a = int(r["host_pos_phage_pos"])
        b = int(r["host_pos_phage_neg"])
        c = int(r["host_neg_phage_pos"])
        d = int(r["host_neg_phage_neg"])
        cis.append(or_ci_haldane(a, b, c, d))

    df["log2OR"] = [c.log2_or for c in cis]
    df["log2OR_lo"] = [c.log2_lo for c in cis]
    df["log2OR_hi"] = [c.log2_hi for c in cis]

    # order: most significant & largest effect
    df = df.sort_values(["q_bh_host_ge", "log2OR"], ascending=[True, False]).reset_index(drop=True)

    # plot top N to keep readable
    top_n = min(30, len(df))
    d = df.head(top_n).copy()

    y = np.arange(len(d))
    sig = d["q_bh_host_ge"].fillna(1.0).to_numpy() < 0.05
    # green-purple scheme (non-sig = green, sig = purple)
    colors = np.where(sig, "#AC99D2", "#70CDBE")

    _set_plot_style()
    fig_h = max(4.5, 0.28 * len(d) + 1.8)
    fig, ax = plt.subplots(figsize=(7.2, fig_h))

    x = d["log2OR"].to_numpy(dtype=float)
    xlo = d["log2OR_lo"].to_numpy(dtype=float)
    xhi = d["log2OR_hi"].to_numpy(dtype=float)
    err_left = x - xlo
    err_right = xhi - x
    ax.errorbar(
        x,
        y,
        xerr=np.vstack([err_left, err_right]),
        fmt="none",
        ecolor="0.25",
        elinewidth=1.3,
        capsize=3,
        capthick=1.0,
        zorder=2,
    )
    ax.scatter(x, y, c=colors, s=72, edgecolors="black", linewidths=0.4, zorder=3)

    ax.axvline(0.0, color="0.2", lw=1.0, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(d["KO"].astype(str).tolist(), fontsize=11)
    ax.set_xlabel(r"log2(OR)  (host KO presence -> phage carries same AMG KO)", fontsize=12)
    ax.set_title("Directionality of KO coupling on CRISPR-linked pairs", fontweight="bold", fontsize=13)
    ax.grid(axis="x", linestyle="--", alpha=0.22, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(axis="x", labelsize=11)
    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def fig_y_bar_bridge_kos(path: Path, out_png: Path, out_pdf: Path) -> None:
    df = pd.read_csv(path, sep="\t", dtype={"KO": str})
    if df.empty:
        return

    df["n_host_genera_spanned"] = pd.to_numeric(df["n_host_genera_spanned"], errors="coerce")
    df["n_linked_MAGs"] = pd.to_numeric(df["n_linked_MAGs"], errors="coerce")
    df["n_linked_vOTUs"] = pd.to_numeric(df["n_linked_vOTUs"], errors="coerce")
    df = df.dropna(subset=["n_host_genera_spanned", "n_linked_MAGs", "n_linked_vOTUs"]).copy()

    top_n = 25
    d = df.sort_values(["n_host_genera_spanned", "n_linked_vOTUs"], ascending=[False, False]).head(top_n).copy()
    d = d.sort_values("n_host_genera_spanned", ascending=True).reset_index(drop=True)

    _set_plot_style()
    fig, ax = plt.subplots(figsize=(9.2, max(5.2, 0.32 * len(d) + 2.2)))

    # color by linked vOTUs (log-scale) to reduce range compression
    v = d["n_linked_vOTUs"].to_numpy(dtype=float)
    norm = mpl.colors.Normalize(vmin=np.log10(max(1.0, v.min())), vmax=np.log10(max(1.0, v.max())))
    # unified green-purple scheme
    cmap = LinearSegmentedColormap.from_list("bridge", ["#70CDBE", "#AC99D2"])
    cols = [cmap(norm(math.log10(max(1.0, x)))) for x in v]

    y = np.arange(len(d))
    ax.barh(y, d["n_host_genera_spanned"].to_numpy(dtype=float), color=cols, edgecolor="white", linewidth=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(d["KO"].astype(str).tolist(), fontsize=9)
    ax.set_xlabel("Number of host genera spanned (CRISPR-linked)")
    ax.set_title("Top bridge KOs across host genera", fontweight="bold")

    # annotate MAG/vOTU coverage
    for i, r in d.iterrows():
        ax.text(
            float(r["n_host_genera_spanned"]) + 2.0,
            i,
            f"MAG={int(r['n_linked_MAGs'])}, vOTU={int(r['n_linked_vOTUs'])}",
            va="center",
            fontsize=8,
            color="0.25",
        )

    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("log10(linked vOTUs)")
    ax.grid(axis="x", linestyle="--", alpha=0.22, linewidth=0.8)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def fig_z_bipartite_network(path: Path, out_png: Path, out_pdf: Path) -> None:
    df = pd.read_csv(path, sep="\t", dtype={"KO": str, "host_genera": str})
    if df.empty:
        return

    df["n_host_genera_spanned"] = pd.to_numeric(df["n_host_genera_spanned"], errors="coerce")
    df["n_linked_MAGs"] = pd.to_numeric(df["n_linked_MAGs"], errors="coerce")
    df["n_linked_vOTUs"] = pd.to_numeric(df["n_linked_vOTUs"], errors="coerce")
    df = df.dropna(subset=["n_host_genera_spanned", "n_linked_MAGs", "n_linked_vOTUs"]).copy()

    top_kos = 15
    d = df.sort_values(["n_host_genera_spanned", "n_linked_vOTUs"], ascending=[False, False]).head(top_kos).copy()

    ko_nodes = d["KO"].astype(str).tolist()
    ko_meta = d.set_index("KO")[["n_linked_MAGs", "n_linked_vOTUs", "n_host_genera_spanned"]].to_dict(orient="index")

    edges: list[tuple[str, str]] = []
    genus_degree: dict[str, int] = defaultdict(int)
    for _, r in d.iterrows():
        ko = str(r["KO"])
        gens = [g.strip() for g in str(r["host_genera"]).split(";") if g.strip()]
        for g in gens:
            edges.append((ko, g))
            genus_degree[g] = genus_degree.get(g, 0) + 1

    # keep only top genera to avoid an unreadable hairball
    max_genera = 60
    top_genera = [g for g, _ in sorted(genus_degree.items(), key=lambda x: (-x[1], x[0]))[:max_genera]]
    edges = [(ko, g) for ko, g in edges if g in set(top_genera)]

    # build graph
    B = nx.Graph()
    B.add_nodes_from(ko_nodes, bipartite=0)
    B.add_nodes_from(top_genera, bipartite=1)
    B.add_edges_from(edges)

    # layout: two columns
    ko_pos = {ko: (0.0, i) for i, ko in enumerate(ko_nodes)}
    gen_sorted = sorted(top_genera, key=lambda g: (-genus_degree.get(g, 0), g))
    gen_pos = {g: (1.0, i) for i, g in enumerate(gen_sorted)}
    pos = {**ko_pos, **gen_pos}

    _set_plot_style()
    fig_h = max(7.0, 0.23 * max(len(ko_nodes), len(gen_sorted)) + 3.0)
    fig, ax = plt.subplots(figsize=(12.5, fig_h))

    # node sizes
    ko_sizes = []
    for ko in ko_nodes:
        m = ko_meta.get(ko, {})
        s = 120 + 0.06 * float(m.get("n_linked_vOTUs", 0)) + 0.03 * float(m.get("n_linked_MAGs", 0))
        ko_sizes.append(s)
    gen_sizes = [60 + 70 * genus_degree.get(g, 0) for g in gen_sorted]

    nx.draw_networkx_edges(B, pos, ax=ax, width=0.65, alpha=0.22, edge_color="0.25")
    nx.draw_networkx_nodes(
        B,
        pos,
        nodelist=ko_nodes,
        node_size=ko_sizes,
        node_color=PALETTE_HEX[5],
        edgecolors="black",
        linewidths=0.45,
        ax=ax,
    )
    nx.draw_networkx_nodes(
        B,
        pos,
        nodelist=gen_sorted,
        node_size=gen_sizes,
        node_color=PALETTE_HEX[0],
        edgecolors="black",
        linewidths=0.35,
        ax=ax,
    )

    nx.draw_networkx_labels(B, pos, labels={k: k for k in ko_nodes}, font_size=9, ax=ax)
    # genus labels smaller
    nx.draw_networkx_labels(B, pos, labels={g: g for g in gen_sorted}, font_size=7, ax=ax)

    ax.set_title("Bridge KO–host genus bipartite network (top bridge KOs)", fontweight="bold")
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def _circular_mean_angle(thetas: np.ndarray) -> float:
    if thetas.size == 0:
        return 0.0
    return float(np.angle(np.mean(np.exp(1j * thetas.astype(float)))))


def fig_coupling_breadth_bipartite(
    bridge_path: Path,
    fisher_path: Path,
    ko_genus_counts_path: Path,
    out_png: Path,
    out_pdf: Path,
    *,
    top_kos: int = 12,
    max_genera: int = 55,
    edge_lw_min: float = 0.42,
    edge_lw_max: float = 1.55,
    edge_alpha_min: float = 0.14,
    edge_alpha_max: float = 0.44,
) -> None:
    """
    Bipartite network on concentric circles (cool-toned palette):
    - Inner ring: AMG KOs (cyan fill, Fisher-colored rim); node size ~ taxonomic breadth.
    - Outer ring: host genera (blue–slate); ordered by mean angle of linked KOs to reduce edge crossing.
    - Edges: chords; width in [edge_lw_min, edge_lw_max], alpha in [edge_alpha_min, edge_alpha_max],
      scaled by vOTU support when `vibrant_amg_KO_by_host_genus.tsv` exists.
    No extra inputs beyond the existing bridge / Fisher / optional KO–genus count tables.
    """
    df = pd.read_csv(bridge_path, sep="\t", dtype={"KO": str, "host_genera": str})
    if df.empty:
        return
    df["n_host_genera_spanned"] = pd.to_numeric(df["n_host_genera_spanned"], errors="coerce")
    df["n_linked_MAGs"] = pd.to_numeric(df["n_linked_MAGs"], errors="coerce")
    df["n_linked_vOTUs"] = pd.to_numeric(df["n_linked_vOTUs"], errors="coerce")
    df = df.dropna(subset=["n_host_genera_spanned", "n_linked_MAGs", "n_linked_vOTUs"]).copy()

    d = df.sort_values(["n_host_genera_spanned", "n_linked_vOTUs"], ascending=[False, False]).head(top_kos).copy()
    ko_nodes = d["KO"].astype(str).tolist()
    ko_meta = d.set_index("KO")[
        ["n_linked_MAGs", "n_linked_vOTUs", "n_host_genera_spanned"]
    ].to_dict(orient="index")

    fisher_map: dict[str, dict[str, float]] = {}
    if fisher_path.exists():
        df_f = pd.read_csv(fisher_path, sep="\t", dtype={"KO": str})
        for _, r in df_f.iterrows():
            ko = str(r["KO"])
            fisher_map[ko] = {
                "q": float(r["q_bh_host_ge"]) if pd.notna(r.get("q_bh_host_ge")) else 1.0,
                "or": float(r["odds_ratio_host_ge"]) if pd.notna(r.get("odds_ratio_host_ge")) else 1.0,
            }

    w_map: dict[tuple[str, str], float] = {}
    if ko_genus_counts_path.exists():
        df_vk = pd.read_csv(ko_genus_counts_path, sep="\t", dtype=str)
        if {"host_genus", "AMG_KO", "n_vOTUs"}.issubset(df_vk.columns):
            for _, r in df_vk.iterrows():
                ko = str(r["AMG_KO"]).strip()
                g = str(r["host_genus"]).strip()
                try:
                    nv = float(r["n_vOTUs"])
                except (TypeError, ValueError):
                    continue
                w_map[(ko, g)] = max(w_map.get((ko, g), 0.0), nv)

    edges: list[tuple[str, str]] = []
    genus_degree: dict[str, int] = defaultdict(int)
    for _, r in d.iterrows():
        ko = str(r["KO"])
        gens = [g.strip() for g in str(r["host_genera"]).split(";") if g.strip()]
        for g in gens:
            edges.append((ko, g))
            genus_degree[g] = genus_degree.get(g, 0) + 1

    top_genera = [g for g, _ in sorted(genus_degree.items(), key=lambda x: (-x[1], x[0]))[:max_genera]]
    top_set = set(top_genera)
    edgelist = [(ko, g) for ko, g in edges if g in top_set]

    raw_w = [float(w_map.get((ko, g), 1.0)) for ko, g in edgelist]
    lo, hi = (min(raw_w), max(raw_w)) if raw_w else (1.0, 1.0)
    span = max(hi - lo, 1e-6)
    w_norm = [(w - lo) / span for w in raw_w]

    gen_sorted = sorted(top_genera, key=lambda g: (-genus_degree.get(g, 0), g))

    gen_to_kos: dict[str, list[str]] = defaultdict(list)
    for ko, g in edgelist:
        gen_to_kos[g].append(ko)

    n_ko = len(ko_nodes)
    n_gen = len(gen_sorted)
    rot = math.pi / 2.0
    theta_ko = (np.linspace(0.0, 2.0 * math.pi, n_ko, endpoint=False) + rot).astype(float)
    ko_theta = {ko: float(theta_ko[i]) for i, ko in enumerate(ko_nodes)}

    def _genus_angle_key(g: str) -> float:
        tt = np.array([ko_theta[k] for k in gen_to_kos.get(g, [])], dtype=float)
        return _circular_mean_angle(tt)

    gen_ordered = sorted(gen_sorted, key=_genus_angle_key)
    theta_gen = (np.linspace(0.0, 2.0 * math.pi, n_gen, endpoint=False) + rot).astype(float)

    r_ko, r_gen = 1.15, 2.75
    pos: dict[str, tuple[float, float]] = {}
    for i, ko in enumerate(ko_nodes):
        th = float(theta_ko[i])
        pos[ko] = (r_ko * math.cos(th), r_ko * math.sin(th))
    for j, g in enumerate(gen_ordered):
        th = float(theta_gen[j])
        pos[g] = (r_gen * math.cos(th), r_gen * math.sin(th))

    def _ko_ring(ko: str) -> tuple[str, float]:
        m = fisher_map.get(ko)
        if not m:
            return "#B4B8BC", 1.35
        q, orv = m["q"], m["or"]
        if q < 0.05 and orv > 1.0:
            return "#5A3D7A", 2.05
        if q < 0.05 and orv <= 1.0:
            return "#2A6056", 1.85
        return "#6BA88E", 1.55

    ng_max = max(float(ko_meta[k]["n_host_genera_spanned"]) for k in ko_nodes)
    ko_sizes = np.array(
        [
            320.0 + 900.0 * math.sqrt(float(ko_meta[ko]["n_host_genera_spanned"]) / max(ng_max, 1.0))
            for ko in ko_nodes
        ],
        dtype=float,
    )
    deg_max = max(genus_degree.values()) if genus_degree else 1
    gen_sizes = np.array(
        [180.0 + 520.0 * (genus_degree[g] / max(deg_max, 1)) ** 0.85 for g in gen_ordered],
        dtype=float,
    )

    # --- draw (cool palette, no label outline) ---
    mpl.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.dpi": 320,
            "font.size": 11,
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "Microsoft YaHei"],
            "axes.unicode_minus": False,
            "pdf.fonttype": 42,
        }
    )
    fig_w = 10.8
    bg = "#F6F8FA"
    fig, ax = plt.subplots(figsize=(fig_w, fig_w), facecolor=bg)
    ax.set_facecolor(bg)

    # faint disk so the plot reads as one object on the page
    ax.add_patch(
        Circle(
            (0.0, 0.0),
            r_gen + 0.55,
            facecolor="#FFFFFF",
            edgecolor="#E4E9EE",
            linewidth=0.9,
            zorder=0,
            alpha=1.0,
        )
    )

    for r_ring in (r_gen, r_ko):
        ax.add_patch(
            Circle(
                (0.0, 0.0),
                r_ring,
                fill=False,
                edgecolor="#DDE3EA",
                linewidth=0.75,
                linestyle=(0, (1, 3)),
                alpha=0.85,
                zorder=0,
            )
        )

    lw_lo, lw_hi = (edge_lw_min, edge_lw_max) if edge_lw_max >= edge_lw_min else (edge_lw_max, edge_lw_min)
    al_lo, al_hi = (edge_alpha_min, edge_alpha_max) if edge_alpha_max >= edge_alpha_min else (
        edge_alpha_max,
        edge_alpha_min,
    )
    segments = [[pos[u], pos[v]] for u, v in edgelist]
    edge_lw = np.array([lw_lo + (lw_hi - lw_lo) * wn for wn in w_norm], dtype=float)
    edge_alpha = np.array([al_lo + (al_hi - al_lo) * wn for wn in w_norm], dtype=float)
    edge_colors = [(0.24, 0.32, 0.44, float(a)) for a in edge_alpha]
    if segments:
        lc = LineCollection(
            segments,
            colors=edge_colors,
            linewidths=edge_lw,
            capstyle="round",
            joinstyle="round",
            zorder=1,
        )
        ax.add_collection(lc)

    xko_arr = np.array([pos[k][0] for k in ko_nodes])
    yko_arr = np.array([pos[k][1] for k in ko_nodes])
    ko_face = "#B8DCF0"
    ko_halo = "#EEF6FB"
    ko_ec = [_ko_ring(ko)[0] for ko in ko_nodes]
    ko_lw = [_ko_ring(ko)[1] for ko in ko_nodes]
    ax.scatter(
        xko_arr,
        yko_arr,
        s=ko_sizes * 1.12,
        c=ko_halo,
        edgecolors="#E8F4FA",
        linewidths=0.9,
        zorder=3,
        alpha=1.0,
    )
    ax.scatter(
        xko_arr,
        yko_arr,
        s=ko_sizes,
        c=ko_face,
        edgecolors=ko_ec,
        linewidths=ko_lw,
        zorder=4,
        alpha=1.0,
    )

    xg_arr = np.array([pos[g][0] for g in gen_ordered])
    yg_arr = np.array([pos[g][1] for g in gen_ordered])
    # Outer ring: cool blue–slate (distinct from inner cyan)
    gen_face = "#94A9C9"
    gen_halo = "#E9EEF6"
    gen_edge = "#4A6285"
    ax.scatter(
        xg_arr,
        yg_arr,
        s=gen_sizes * 1.1,
        c=gen_halo,
        edgecolors="#D4DCE8",
        linewidths=0.85,
        zorder=3,
        alpha=1.0,
    )
    ax.scatter(
        xg_arr,
        yg_arr,
        s=gen_sizes,
        c=gen_face,
        edgecolors=gen_edge,
        linewidths=1.0,
        zorder=4,
        alpha=0.98,
    )

    for ko in ko_nodes:
        x, y = pos[ko]
        dist = max(1e-6, math.hypot(x, y))
        ux, uy = x / dist, y / dist
        ax.annotate(
            ko,
            (x, y),
            xytext=(-ux * 16, -uy * 16),
            textcoords="offset points",
            ha="center",
            va="center",
            fontsize=10.2,
            fontweight="semibold",
            color="#0f172a",
            zorder=5,
        )
    _degs = sorted(genus_degree.values(), reverse=True) if genus_degree else [0]
    deg_hi = _degs[11] if len(_degs) > 11 else _degs[-1]
    # Radial outward offset for genus labels (data coords)
    pad_r = 0.13 * r_gen
    pad_r_hub = 0.15 * r_gen

    def _genus_label_rotation(px: float, py: float) -> float:
        """Degrees: text axis along radius toward origin; flip if upside-down."""
        rot = math.degrees(math.atan2(py, px)) + 180.0
        if rot > 180.0:
            rot -= 360.0
        if rot <= -180.0:
            rot += 360.0
        if rot > 90.0:
            rot -= 180.0
        elif rot < -90.0:
            rot += 180.0
        return rot

    for g in gen_ordered:
        x, y = pos[g]
        dist = max(1e-6, math.hypot(x, y))
        ux, uy = x / dist, y / dist
        fs = 9.2 if genus_degree.get(g, 0) >= max(deg_hi, 3) else 8.2
        pr = pad_r_hub if genus_degree.get(g, 0) >= max(deg_hi, 3) else pad_r
        tx, ty = x + ux * pr, y + uy * pr
        rot = _genus_label_rotation(x, y)
        ax.text(
            tx,
            ty,
            g,
            ha="center",
            va="center",
            rotation=rot,
            rotation_mode="anchor",
            fontsize=fs,
            color="#1e293b",
            style="italic",
            zorder=5,
            clip_on=False,
        )

    lim = r_gen + 1.58
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")

    ax.set_title(
        "CRISPR-linked AMG KOs",
        fontsize=15.5,
        fontweight="bold",
        color="#0f172a",
        pad=10,
    )
    ax.text(
        0.5,
        1.02,
        "Circular layout: inner = AMG KO (phage); outer = host genus",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=11.2,
        color="#64748b",
        fontweight="normal",
    )
    ax.text(
        0.5,
        -0.02,
        "KO size ~ genera spanned; genus size ~ links to top KOs.  KO rim = Fisher.  Edges ~ vOTU support.",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=9.4,
        color="#64748b",
    )

    from matplotlib.patches import Patch

    leg_ko = [
        Patch(facecolor=ko_face, edgecolor="#5A3D7A", linewidth=2.0, label="FDR<0.05, OR>1"),
        Patch(facecolor=ko_face, edgecolor="#2A6056", linewidth=1.85, label="FDR<0.05, OR≤1"),
        Patch(facecolor=ko_face, edgecolor="#6BA88E", linewidth=1.5, label="FDR≥0.05"),
        Patch(facecolor=ko_face, edgecolor="#B4B8BC", linewidth=1.35, label="Not in Fisher table"),
    ]
    leg_host = [Patch(facecolor=gen_face, edgecolor=gen_edge, linewidth=1.0, label="Host genus (hub)")]
    leg = leg_ko + leg_host
    ax.legend(
        handles=leg,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.03),
        bbox_transform=ax.transAxes,
        ncol=3,
        fontsize=9.2,
        frameon=True,
        fancybox=True,
        framealpha=0.96,
        edgecolor="#D8DEE6",
        facecolor="#FFFFFF",
        title="Encoding",
        title_fontsize=9.2,
    )

    fig.tight_layout()
    fig.subplots_adjust(top=0.9, bottom=0.12)
    fig.savefig(out_png, bbox_inches="tight", facecolor=fig.get_facecolor())
    fig.savefig(out_pdf, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    _set_plot_style()

    fig_x_forest_directionality(
        IN_DIR / "ko_directionality_fisher_linked_pairs.tsv",
        OUT_DIR / "FigX_directionality_forest_log2OR.png",
        OUT_DIR / "FigX_directionality_forest_log2OR.pdf",
    )
    fig_y_bar_bridge_kos(
        IN_DIR / "bridge_kos_all_genera_counts.tsv",
        OUT_DIR / "FigY_top_bridge_KOs_bar.png",
        OUT_DIR / "FigY_top_bridge_KOs_bar.pdf",
    )
    fig_z_bipartite_network(
        IN_DIR / "bridge_kos_all_genera_counts.tsv",
        OUT_DIR / "FigZ_bridge_KO_genus_bipartite_network.png",
        OUT_DIR / "FigZ_bridge_KO_genus_bipartite_network.pdf",
    )
    fig_coupling_breadth_bipartite(
        IN_DIR / "bridge_kos_all_genera_counts.tsv",
        IN_DIR / "ko_directionality_fisher_linked_pairs.tsv",
        ROOT / "summary_tables" / "vibrant_amg_spacer_host" / "vibrant_amg_KO_by_host_genus.tsv",
        OUT_DIR / "Fig_coupling_breadth_KO_genus_network.png",
        OUT_DIR / "Fig_coupling_breadth_KO_genus_network.pdf",
    )

    # Additional visuals requested:
    fig_z2_chord_ko_taxonomy(
        ko_table=IN_DIR / "bridge_kos_all_genera_counts.tsv",
        out_png=OUT_DIR / "FigZ2_chord_KO_phylum.png",
        out_pdf=OUT_DIR / "FigZ2_chord_KO_phylum.pdf",
    )
    fig_z3_taxonomy_rings(
        ko_table=IN_DIR / "bridge_kos_all_genera_counts.tsv",
        out_png=OUT_DIR / "FigZ3_taxonomy_rings_KO_coverage.png",
        out_pdf=OUT_DIR / "FigZ3_taxonomy_rings_KO_coverage.pdf",
    )


def _parse_taxonomy_table(path: Path) -> dict[str, dict[str, str]]:
    """
    Return MAG -> ranks mapping with keys: phylum,class,order,family,genus.
    """
    import re

    def _get(tag: str, s: str) -> str:
        m = re.search(rf"{tag}__([^;]+)", s)
        if not m:
            return "Unclassified"
        v = m.group(1).strip()
        return v if v else "Unclassified"

    out: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n\r")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t", 1)
            if len(parts) < 2:
                continue
            mag, tax = parts[0].strip(), parts[1].strip()
            if not mag:
                continue
            out[mag] = {
                "phylum": _get("p", tax),
                "class": _get("c", tax),
                "order": _get("o", tax),
                "family": _get("f", tax),
                "genus": _get("g", tax),
            }
    return out


def _load_phage_votu_kos(path: Path) -> dict[str, set[str]]:
    df = pd.read_csv(path, sep="\t", dtype=str)
    df["votu"] = df["scaffold"].astype(str).str.replace(r"\|provirus_.*$", "", regex=True).str.strip()
    out: dict[str, set[str]] = defaultdict(set)
    for _, r in df.iterrows():
        v = str(r.get("votu", "")).strip()
        if not v:
            continue
        ko = str(r.get("AMG KO", "")).strip().upper()
        ko = ko.removeprefix("KO:").removeprefix("KO").strip()
        if ko.startswith("K") and len(ko) >= 6:
            # Keep only Kxxxxx / Kxxxxxx
            ko = ko.split()[0].split(",")[0]
            out[v].add(ko)
    return dict(out)


def _load_crispr_pairs(path: Path) -> list[tuple[str, str]]:
    df = pd.read_csv(path, sep="\t", dtype=str)
    df["votu"] = df["viral_contig"].astype(str).str.replace(r"\|provirus_.*$", "", regex=True).str.strip()
    df["mag_id"] = df["mag_id"].astype(str).str.strip()
    sub = df.dropna(subset=["mag_id", "votu"]).drop_duplicates(subset=["mag_id", "votu"])
    return [(str(a), str(b)) for a, b in zip(sub["mag_id"], sub["votu"])]


def _compute_ko_genus_counts(top_kos: list[str]) -> tuple[pd.DataFrame, dict[str, dict[str, str]]]:
    """
    Edge weights for selected KOs:
    count = number of unique CRISPR (mag,votu) pairs where votu carries KO and mag maps to genus/phylum.
    """
    tax = _parse_taxonomy_table(TAX_PATH)
    votu_kos = _load_phage_votu_kos(VIB_INDIV)
    pairs = _load_crispr_pairs(CRISPR_BEST)

    top = set(top_kos)
    ko_genus: dict[tuple[str, str], int] = defaultdict(int)
    for mag, votu in pairs:
        ranks = tax.get(mag)
        if not ranks:
            continue
        genus = ranks.get("genus", "Unclassified")
        kos = votu_kos.get(votu)
        if not kos:
            continue
        for ko in kos:
            if ko in top:
                ko_genus[(ko, genus)] += 1

    rows = [{"KO": k, "genus": g, "n_pairs": n} for (k, g), n in ko_genus.items()]
    df = pd.DataFrame(rows)
    return df, tax


def fig_z2_chord_ko_taxonomy(ko_table: Path, out_png: Path, out_pdf: Path) -> None:
    """
    Chord-like diagram: KO ↔ host phylum, weighted by CRISPR-linked pair counts.
    """
    dfk = pd.read_csv(ko_table, sep="\t", dtype=str)
    if dfk.empty:
        return
    dfk["n_host_genera_spanned"] = pd.to_numeric(dfk["n_host_genera_spanned"], errors="coerce")
    dfk = dfk.dropna(subset=["n_host_genera_spanned"]).sort_values("n_host_genera_spanned", ascending=False)
    top_kos = dfk["KO"].astype(str).head(12).tolist()

    df_edges, tax = _compute_ko_genus_counts(top_kos)
    if df_edges.empty:
        return

    # map genus->phylum using taxonomy table (MAG-based), best-effort via majority vote across MAGs
    genus_to_phylum: dict[str, str] = {}
    tmp: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for mag, ranks in tax.items():
        g = ranks.get("genus", "Unclassified")
        p = ranks.get("phylum", "Unclassified")
        tmp[g][p] += 1
    for g, c in tmp.items():
        genus_to_phylum[g] = max(c.items(), key=lambda x: x[1])[0]

    df_edges["phylum"] = df_edges["genus"].map(lambda g: genus_to_phylum.get(str(g), "Unclassified"))
    df_pw = df_edges.groupby(["KO", "phylum"], as_index=False)["n_pairs"].sum()

    # pick top phyla to keep readable
    phylum_order = (
        df_pw.groupby("phylum")["n_pairs"].sum().sort_values(ascending=False).head(12).index.astype(str).tolist()
    )
    df_pw = df_pw[df_pw["phylum"].isin(phylum_order)].copy()

    # node order
    kos = top_kos
    phyla = phylum_order
    w = (
        df_pw.pivot_table(index="KO", columns="phylum", values="n_pairs", fill_value=0, aggfunc="sum")
        .reindex(index=kos, columns=phyla)
        .to_numpy(dtype=float)
    )
    if w.sum() <= 0:
        return

    # angles for KO arc (left) and phylum arc (right)
    _set_plot_style()
    fig = plt.figure(figsize=(10.5, 9.5))
    ax = fig.add_subplot(111, polar=True)
    ax.set_axis_off()

    def _arc_positions(n: int, start: float, end: float) -> list[float]:
        if n == 1:
            return [(start + end) / 2]
        return list(np.linspace(start, end, n))

    ko_angles = _arc_positions(len(kos), math.radians(110), math.radians(250))
    ph_angles = _arc_positions(len(phyla), math.radians(-70), math.radians(70))

    # draw labels and small wedges
    r0 = 1.0
    for ang, label in zip(ko_angles, kos):
        ax.text(ang, r0 + 0.08, label, ha="center", va="center", fontsize=9, rotation=np.degrees(ang) - 180)
    for ang, label in zip(ph_angles, phyla):
        ax.text(ang, r0 + 0.10, label, ha="center", va="center", fontsize=9, rotation=np.degrees(ang))

    # chords as Bezier curves in polar coordinates projected to cartesian
    vmax = w.max()
    for i, ko in enumerate(kos):
        for j, ph in enumerate(phyla):
            ww = w[i, j]
            if ww <= 0:
                continue
            a1, a2 = ko_angles[i], ph_angles[j]
            # convert to cartesian coordinates
            x1, y1 = r0 * math.cos(a1), r0 * math.sin(a1)
            x2, y2 = r0 * math.cos(a2), r0 * math.sin(a2)
            # control points pull toward center
            c = 0.15
            cx1, cy1 = c * math.cos(a1), c * math.sin(a1)
            cx2, cy2 = c * math.cos(a2), c * math.sin(a2)
            path = MplPath(
                [(x1, y1), (cx1, cy1), (cx2, cy2), (x2, y2)],
                [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4],
            )
            alpha = 0.12 + 0.65 * (ww / vmax)
            lw = 0.4 + 2.2 * (ww / vmax)
            patch = PathPatch(path, facecolor="none", edgecolor="#2c7bb6", lw=lw, alpha=alpha)
            ax.add_patch(patch)

    ax.set_title("Chord-like KO <-> host phylum links (weighted by CRISPR-linked pairs)", pad=28, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def fig_z3_taxonomy_rings(ko_table: Path, out_png: Path, out_pdf: Path) -> None:
    """
    Circular taxonomy-ordered 'tree' + KO coverage rings.
    Note: this is taxonomy-ordered (not a true phylogenetic tree).
    """
    dfk = pd.read_csv(ko_table, sep="\t", dtype=str)
    if dfk.empty:
        return
    dfk["n_host_genera_spanned"] = pd.to_numeric(dfk["n_host_genera_spanned"], errors="coerce")
    dfk = dfk.dropna(subset=["n_host_genera_spanned"]).sort_values("n_host_genera_spanned", ascending=False)
    top_kos = dfk["KO"].astype(str).head(6).tolist()

    df_edges, tax = _compute_ko_genus_counts(top_kos)
    if df_edges.empty:
        return

    # Build genus -> phylum mapping (majority vote)
    genus_to_phylum: dict[str, str] = {}
    tmp: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for mag, ranks in tax.items():
        g = ranks.get("genus", "Unclassified")
        p = ranks.get("phylum", "Unclassified")
        tmp[g][p] += 1
    for g, c in tmp.items():
        genus_to_phylum[g] = max(c.items(), key=lambda x: x[1])[0]

    # select genera that appear in edges
    genera = sorted(
        df_edges["genus"].unique().tolist(),
        key=lambda g: (genus_to_phylum.get(g, "Unclassified"), g),
    )
    if len(genera) < 10:
        return

    # matrix: genus x KO = n_pairs
    mat = (
        df_edges.pivot_table(index="genus", columns="KO", values="n_pairs", fill_value=0, aggfunc="sum")
        .reindex(index=genera, columns=top_kos)
        .to_numpy(dtype=float)
    )
    mat = np.log1p(mat)  # compress

    # phylum segments
    phylum_of = [genus_to_phylum.get(g, "Unclassified") for g in genera]
    phylum_blocks: list[tuple[str, int, int]] = []
    start = 0
    for i in range(1, len(genera) + 1):
        if i == len(genera) or phylum_of[i] != phylum_of[start]:
            phylum_blocks.append((phylum_of[start], start, i))
            start = i

    _set_plot_style()
    fig = plt.figure(figsize=(12.5, 12.5))
    ax = fig.add_subplot(111, polar=True)
    ax.set_axis_off()

    n = len(genera)
    angles = np.linspace(0, 2 * math.pi, n, endpoint=False)
    width = 2 * math.pi / n

    # inner taxonomy ring (phylum background)
    phyla = [b[0] for b in phylum_blocks]
    ph_cmap = mpl.colormaps.get_cmap("tab20")
    ph_color = {p: ph_cmap(i % 20) for i, p in enumerate(phyla)}
    r_inner = 0.55
    ring_thick = 0.10
    for p, a, b in phylum_blocks:
        ang0 = angles[a]
        ang1 = angles[b - 1] + width
        wedge = Wedge((0, 0), r_inner + ring_thick, np.degrees(ang0), np.degrees(ang1), width=ring_thick)
        wedge.set_facecolor(ph_color.get(p))
        wedge.set_edgecolor("white")
        wedge.set_linewidth(0.3)
        ax.add_patch(wedge)

    # KO rings
    cmap = mpl.colormaps.get_cmap("YlOrRd")
    vmin, vmax = float(np.min(mat)), float(np.max(mat))
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax if vmax > vmin else (vmin + 1.0))

    base = 0.70
    step = 0.07
    for ki, ko in enumerate(top_kos):
        r0 = base + ki * step
        for gi in range(n):
            val = mat[gi, ki]
            col = cmap(norm(val))
            wedge = Wedge((0, 0), r0 + step * 0.90, np.degrees(angles[gi]), np.degrees(angles[gi] + width), width=step * 0.90)
            wedge.set_facecolor(col)
            wedge.set_edgecolor("none")
            ax.add_patch(wedge)
        ax.text(math.radians(90), r0 + step * 0.45, ko, ha="center", va="center", fontsize=9)

    # genus tick labels (sparse)
    label_every = max(1, n // 40)
    for i, (ang, g) in enumerate(zip(angles, genera)):
        if i % label_every != 0:
            continue
        ax.text(ang, base + len(top_kos) * step + 0.06, g, fontsize=6, ha="center", va="center", rotation=np.degrees(ang))

    ax.set_title(
        "Fig Z3. Taxonomy-ordered circular map + KO coverage rings\n"
        "Inner: host phylum blocks; outer rings: log1p(#CRISPR-linked pairs carrying KO)",
        pad=30,
    )
    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()

