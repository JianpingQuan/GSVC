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


"""
Python: beta diversity — geographic grouping (country or continent), centroid position + within-group dispersion plots + PERMANOVA / betadisper

Visualization outputs (see end of main):
  - PCoA (Aitchison / Bray) colored by `--geo` level + **group centroids**
  - Extra: PCoA **country=color, continent=marker** (dual legend; lists **all countries** by default; use `--country-plot-top K` to keep top K countries + Other)
  - **Within-group dispersion**: distance from samples to group centroid in PCoA space (first r dimensions) — **boxplot**

Dependencies: numpy, pandas, scipy, scikit-learn, matplotlib

Run examples:
  python compositional_continent_analysis.py --geo continent
  python compositional_continent_analysis.py --geo country
  python compositional_continent_analysis.py --geo continent --max-taxa 15000 --max-samples 800

Environment variable PHAGE_PROJECT_ROOT: Phage project root directory

Server batch-run instructions: see SERVER_RUNBOOK.md in this directory
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.manifold import MDS
from sklearn.metrics import pairwise_distances

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
_THIS = Path(__file__).resolve()
_DEFAULT_ROOT = _THIS.parent.parent
ROOT = Path(os.environ.get("PHAGE_PROJECT_ROOT", str(_DEFAULT_ROOT))).expanduser().resolve()
if not ROOT.exists():
    ROOT = _DEFAULT_ROOT

sys.path.insert(0, str(ROOT))
from biogeography_descriptive_analysis import (  # noqa: E402
    MATRIX,
    continent_from_country,
    load_meta,
)

try:
    from biogeography_descriptive_analysis import (  # noqa: E402
        apply_matplotlib_pdf_editable_arial,
    )
except ImportError:
    # If server lacks updated biogeography_descriptive_analysis.py, use local equivalent
    def apply_matplotlib_pdf_editable_arial(base_size: float = 10.5) -> None:
        import matplotlib as mpl

        tick = max(8.0, base_size - 1.0)
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
                "font.size": base_size,
                "axes.titlesize": base_size + 1.0,
                "axes.labelsize": base_size,
                "xtick.labelsize": tick,
                "ytick.labelsize": tick,
                "legend.fontsize": tick,
            }
        )

OUT_DIR = ROOT / "04_CLR_PERMANOVA"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PREVALENCE_MIN = 0.01
PSEUDOCOUNT = 1.0
PERM = 999
MIN_SAMPLE_SUM = 1
RNG_SEED = 42


def clr_rows(Y: np.ndarray, pseudocount: float = 1.0) -> np.ndarray:
    Z = Y.astype(np.float64) + pseudocount
    logz = np.log(Z)
    return logz - logz.mean(axis=1, keepdims=True)


def euclidean_dist_sq_fast(X: np.ndarray) -> np.ndarray:
    """Samples as rows; D_ij = ||x_i - x_j|| (O(n^2 p) but BLAS matrix multiply, much faster than pairwise pdist)."""
    X = np.asarray(X, dtype=np.float64)
    g = X @ X.T
    s = np.sum(X * X, axis=1, keepdims=True)
    d2 = s + s.T - 2.0 * g
    np.clip(d2, 0.0, None, out=d2)
    return np.sqrt(d2, dtype=np.float64)


def cmdscale_pcoa(D: np.ndarray, k: int = 10):
    """Classical MDS / PCoA from distance matrix D (n x n)."""
    D = np.asarray(D, dtype=np.float64)
    n = D.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (D**2) @ J
    w, V = np.linalg.eigh(B)
    idx = np.argsort(w)[::-1]
    w, V = w[idx], V[:, idx]
    w_pos = np.maximum(w[:k], 0.0)
    X = V[:, :k] * np.sqrt(w_pos)
    var_explained = w_pos / max(w_pos.sum(), 1e-12) * 100.0
    return X, var_explained




def permanova_pseudof(D: np.ndarray, group_codes: np.ndarray) -> float:
    """Anderson 2001 pseudo-F; D square symmetric; group_codes int 0..G-1."""
    D = np.asarray(D, dtype=np.float64)
    n = D.shape[0]
    D2 = D**2
    s_T = D2.sum() / n / 2.0
    G = int(group_codes.max()) + 1
    s_W = 0.0
    for g in range(G):
        idx = np.where(group_codes == g)[0]
        ng = len(idx)
        if ng < 2:
            continue
        sub = D2[np.ix_(idx, idx)]
        s_W += sub.sum() / ng / 2.0
    s_A = s_T - s_W
    if G < 2 or n <= G:
        return float("nan")
    return (s_A / (G - 1)) / (s_W / (n - G))


def stratified_shuffle_str(labels: np.ndarray, strata: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    out = np.array(labels, dtype=object, copy=True)
    for s in np.unique(strata):
        m = strata == s
        block = out[m]
        if len(block) <= 1:
            continue
        out[m] = rng.permutation(block)
    return out


def permanova_stratified(
    D: np.ndarray,
    cont: np.ndarray,
    study: np.ndarray,
    categories: list | pd.Index,
    perm: int,
    rng: np.random.Generator,
) -> tuple[float, float, np.ndarray]:
    """Shuffle continent labels within Study; fixed category levels for consistent pseudo-F degrees of freedom."""
    codes = pd.Categorical(cont, categories=categories).codes
    obs = permanova_pseudof(D, codes)
    stats = np.empty(perm, dtype=np.float64)
    for i in range(perm):
        sh = stratified_shuffle_str(cont, study, rng)
        c2 = pd.Categorical(sh, categories=categories).codes
        stats[i] = permanova_pseudof(D, c2)
    p = (1 + np.sum(stats >= obs)) / (1 + perm)
    return obs, float(p), stats


def permanova_unstratified(
    D: np.ndarray, labels: np.ndarray, categories: list | pd.Index, perm: int, rng: np.random.Generator
) -> tuple[float, float]:
    codes = pd.Categorical(labels, categories=categories).codes
    obs = permanova_pseudof(D, codes)
    stats = np.empty(perm, dtype=np.float64)
    for i in range(perm):
        sh = rng.permutation(labels)
        c2 = pd.Categorical(sh, categories=categories).codes
        stats[i] = permanova_pseudof(D, c2)
    p = (1 + np.sum(stats >= obs)) / (1 + perm)
    return obs, float(p)


def _dispersion_to_centroid(X: np.ndarray, group: np.ndarray) -> np.ndarray:
    n = X.shape[0]
    dist_to_centroid = np.zeros(n, dtype=np.float64)
    for g in np.unique(group):
        m = group == g
        if np.sum(m) < 2:
            dist_to_centroid[m] = 0.0
            continue
        cent = X[m].mean(axis=0)
        dist_to_centroid[m] = np.linalg.norm(X[m] - cent, axis=1)
    return dist_to_centroid


def collapse_country_for_plot(countries: np.ndarray, top_n: int, other_label: str = "Other (rare country)") -> np.ndarray:
    """Visualization only: if top_n<=0, no collapse, list all countries in data; if top_n>0, keep top top_n countries by sample size, merge rest into other_label."""
    s = pd.Series(np.asarray(countries).ravel(), dtype="object").astype(str).str.strip()
    s = s.mask(s.str.lower().isin(["", "nan", "na", "none"]), "NA")
    if top_n is None or int(top_n) <= 0:
        return s.fillna("NA").values
    top_n = int(top_n)
    vc = s.value_counts()
    if len(vc) <= top_n:
        return s.fillna("NA").values
    keep = set(vc.head(top_n).index)
    out = s.where(s.isin(keep), other_label)
    return out.fillna("NA").values


def plot_pcoa_country_color_continent_shape(
    x1: np.ndarray,
    x2: np.ndarray,
    country: np.ndarray,
    continent: np.ndarray,
    title: str,
    out_base: Path,
    country_plot_top: int,
) -> None:
    """PCoA: country -> color; continent -> marker (dual legend)."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    co_plot = collapse_country_for_plot(country, country_plot_top)
    cont = np.asarray(continent).astype(str)
    country_cats = sorted(pd.unique(co_plot))
    cont_cats = sorted(pd.unique(cont))
    markers = ("o", "^", "s", "D", "v", "P", "X", "<", ">", "h", "8", "*", "p", "H")
    cont_m = {c: markers[i % len(markers)] for i, c in enumerate(cont_cats)}
    n_co = len(country_cats)
    base_colors = plt.cm.tab20(np.linspace(0, 1, min(20, max(n_co, 2))))
    if n_co > 20:
        extra = plt.cm.Set3(np.linspace(0, 1, n_co - 20))
        base_colors = np.vstack([base_colors, extra])
    co_color = {c: base_colors[i % len(base_colors)] for i, c in enumerate(country_cats)}

    fig_w = 10.5 if n_co <= 14 else (12.0 if n_co <= 22 else 13.5)
    fig_h = 7.0 if n_co <= 14 else 7.8
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    for ct in cont_cats:
        mk = cont_m[ct]
        for co in country_cats:
            m = (co_plot == co) & (cont == ct)
            if not np.any(m):
                continue
            ax.scatter(
                x1[m],
                x2[m],
                s=22,
                alpha=0.55,
                c=[co_color[co]],
                marker=mk,
                edgecolors="white",
                linewidths=0.25,
                label=None,
            )

    leg_country = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=co_color[c],
            markersize=7,
            linestyle="None",
            label=str(c)[:42],
        )
        for c in country_cats
    ]
    leg_cont = [
        Line2D(
            [0],
            [0],
            marker=cont_m[c],
            color="gray",
            markerfacecolor="lightgray",
            markeredgecolor="gray",
            markersize=8,
            linestyle="None",
            label=str(c)[:32],
        )
        for c in cont_cats
    ]
    leg_ncol = 2 if n_co >= 10 else 1
    leg_fs = 8.0 if n_co <= 12 else 7.0
    leg1 = ax.legend(
        handles=leg_country,
        title="Country",
        bbox_to_anchor=(1.02, 1.0),
        loc="upper left",
        fontsize=leg_fs,
        ncol=leg_ncol,
        frameon=False,
    )
    ax.add_artist(leg1)
    ax.legend(
        handles=leg_cont,
        title="Continent",
        bbox_to_anchor=(1.02, 0.35),
        loc="upper left",
        fontsize=9,
        frameon=False,
    )
    ax.set_xlabel("PCoA axis 1")
    ax.set_ylabel("PCoA axis 2")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    for ext in (".png", ".pdf"):
        fig.savefig(out_base.with_suffix(ext), dpi=300 if ext == ".png" else None, bbox_inches="tight")
    plt.close(fig)


def betadisper_perm(D: np.ndarray, group: np.ndarray, perm: int, rng: np.random.Generator) -> tuple[float, float]:
    """Compute PCoA coordinates once; permutations only recompute distances to centroid and F."""
    from scipy import stats as scipy_stats

    n = D.shape[0]
    r = min(8, n - 1)
    X, _ = cmdscale_pcoa(D, k=r)
    dist_to_centroid = _dispersion_to_centroid(X, group)
    groups_list = [dist_to_centroid[group == g] for g in np.unique(group)]
    if len(groups_list) < 2:
        return float("nan"), float("nan")
    f_obs, _ = scipy_stats.f_oneway(*groups_list)
    stats = np.empty(perm, dtype=np.float64)
    grp = np.asarray(group)
    for i in range(perm):
        gperm = rng.permutation(grp)
        d2 = _dispersion_to_centroid(X, gperm)
        gl = [d2[gperm == g] for g in np.unique(grp)]
        f, _ = scipy_stats.f_oneway(*gl)
        stats[i] = f
    p = (1 + np.sum(stats >= f_obs)) / (1 + perm)
    return float(f_obs), float(p)


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Composition+ PERMANOVA（Python）")
    parser.add_argument(
        "--max-taxa",
        type=int,
        default=None,
        metavar="N",
        help="Read only the first N vOTU rows of the matrix for quick plotting (do not use for formal analysis)",
    )
    parser.add_argument(
        "--perm",
        type=int,
        default=None,
        metavar="N",
        help=f"Number of permutations; default {PERM}; default is 99 when used with --max-taxa.",
    )
    parser.add_argument(
        "--nmds",
        action="store_true",
        help="Force NMDS calculation (very slow for large samples; skipped by default in preview mode)",
    )
    parser.add_argument(
        "--max-samples",
        type=int,
        default=None,
        metavar="M",
        help="Use at most M samples (random subset); defaults to 800 when used with --max-taxa and not specified.",
    )
    parser.add_argument(
        "--geo",
        choices=("continent", "country"),
        default="continent",
        help="Geographic grouping: continent or country; used for coloring in PERMANOVA, betadisper, and primary PCoA.",
    )
    parser.add_argument(
        "--country-plot-top",
        type=int,
        default=0,
        metavar="K",
        help="Dual-coding PCoA: When K=0 (default), all countries are listed; when K>0, only the top K countries (by sample size) are individually colored, while the rest are grouped as "Other" (this does not affect the statistical test).",
    )
    parser.add_argument(
        "--no-dual-pcoa",
        action="store_true",
        help="Do not output Aitchison/Bray PCoA plots using "country color + continent shape" coding.",
    )
    args = parser.parse_args()
    max_taxa = args.max_taxa
    if args.perm is not None:
        n_perm = max(29, int(args.perm))
    elif max_taxa is not None:
        n_perm = 99
    else:
        n_perm = PERM
    max_samples = args.max_samples
    if max_taxa is not None and max_samples is None:
        max_samples = 800
    suffix = f"_previewTaxa{max_taxa}" if max_taxa is not None else ""
    if max_samples is not None:
        suffix = f"{suffix}_previewS{max_samples}"
    suffix_full = f"{suffix}_geo{args.geo}"
    run_nmds = args.nmds or (max_taxa is None)

    rng = np.random.default_rng(RNG_SEED)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt
    from matplotlib.transforms import blended_transform_factory

    apply_matplotlib_pdf_editable_arial(10.5)

    def qualitative_colors_n(n: int) -> np.ndarray:
        """Qualitative group colors: <=10 use tab10; <=20 use tab20 interleaved hues; more classes use uniform HSV hues to avoid adjacent tab20 colors being too similar."""
        n = max(int(n), 1)
        if n <= 10:
            return plt.cm.tab10(np.linspace(0, 1, n, endpoint=False))
        if n <= 20:
            t = plt.cm.tab20(np.linspace(0, 1, 20))
            order = list(range(0, 20, 2)) + list(range(1, 20, 2))
            return t[order[:n]]
        hues = np.linspace(0, 1, n, endpoint=False)
        return np.array([mcolors.hsv_to_rgb((float(h), 0.72, 0.88)) for h in hues])

    print("ROOT =", ROOT)
    if max_taxa is not None or max_samples is not None:
        print(
            "preview: max_taxa=",
            max_taxa,
            " max_samples=",
            max_samples,
            " perm=",
            n_perm,
            " NMDS=",
            run_nmds,
            flush=True,
        )
    meta = load_meta()
    meta = meta.drop_duplicates(subset=["sample_id"], keep="first")
    meta["continent"] = meta["Country"].map(lambda x: continent_from_country(str(x)))
    meta["Study"] = meta["Study"].astype(str)

    # Read matrix header, load only matching columns
    with MATRIX.open("r", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
    if not header or header[0] != "Contig":
        raise ValueError("First column should be Contig")
    count_cols = header[1:]
    strip = lambda c: c[: -len("_Count")] if c.endswith("_Count") else c
    id_from_col = [strip(c) for c in count_cols]
    id_set = set(id_from_col)
    seen: set[str] = set()
    want_ids: list[str] = []
    for s in meta["sample_id"].astype(str):
        if s in id_set and s not in seen:
            seen.add(s)
            want_ids.append(s)
    usecols = ["Contig"] + [count_cols[id_from_col.index(s)] for s in want_ids]

    print("Read the number of matrix columns:", len(usecols) - 1, flush=True)
    dtype_map = {c: np.int32 for c in usecols[1:]}
    read_kw = dict(
        sep="\t",
        usecols=usecols,
        dtype=dtype_map,
        converters={"Contig": str},
        low_memory=False,
    )
    if max_taxa is not None:
        read_kw["nrows"] = int(max_taxa)
    mat = pd.read_csv(MATRIX, **read_kw)
    mat = mat.set_index("Contig")
    mat.columns = want_ids
    sid_order = list(mat.columns)
    meta["sample_id"] = meta["sample_id"].astype(str)
    meta = meta.drop_duplicates(subset=["sample_id"], keep="first")
    meta = meta[meta["sample_id"].isin(sid_order)].set_index("sample_id").loc[sid_order].reset_index()

    Y = mat.values.T.astype(np.float64)
    sample_ids = sid_order
    depths = Y.sum(axis=1)
    keep = depths >= MIN_SAMPLE_SUM
    Y = Y[keep]
    meta = meta.iloc[keep].reset_index(drop=True)
    sample_ids = [sample_ids[i] for i in np.flatnonzero(keep)]

    prev = (Y > 0).mean(axis=0)
    Y = Y[:, prev >= PREVALENCE_MIN]
    print("Filter:", Y.shape[0], "Sample,", Y.shape[1], "taxa", flush=True)

    if max_samples is not None and Y.shape[0] > max_samples:
        sub = rng.choice(Y.shape[0], size=int(max_samples), replace=False)
        sub.sort()
        Y = Y[sub]
        meta = meta.iloc[sub].reset_index(drop=True)
        sample_ids = [sample_ids[i] for i in sub]
        print("After subsampling", Y.shape[0], "Sample", flush=True)

    clr_mat = clr_rows(Y, PSEUDOCOUNT)
    print("Calculate the Aitchison distance matrix …", flush=True)
    D_ait = euclidean_dist_sq_fast(clr_mat)

    row_sum = Y.sum(axis=1, keepdims=True)
    row_sum[row_sum == 0] = 1.0
    Y_rel = Y / row_sum
    print("Calculating the Bray-Curtis distance matrix …", flush=True)
    # Bray on full dimensions is very slow for high dim + large n; preview uses high-abundance column subset
    if max_taxa is not None:
        n_br = min(500, Y_rel.shape[1])
        top_tax = np.argsort(-Y_rel.sum(axis=0))[:n_br]
        Y_br = Y_rel[:, top_tax]
        print("  (Preview: Bray-Curtis distance based on the top ", n_br, " most abundant vOTUs)", flush=True)
    else:
        Y_br = Y_rel
    # n_jobs=-1 can hang on Windows; single-thread is more stable
    D_bray = pairwise_distances(Y_br, metric="braycurtis", n_jobs=1)

    study = meta["Study"].astype(str).values
    study_cats = sorted(pd.unique(study))
    if args.geo == "country":
        geo = meta["Country"].astype(str).str.strip().values
    else:
        geo = meta["continent"].astype(str).values
    geo_cats = sorted(pd.unique(geo))

    np.save(OUT_DIR / f"dist_aitchison{suffix_full}.npy", D_ait)
    np.save(OUT_DIR / f"dist_bray_curtis{suffix_full}.npy", D_bray)
    if max_taxa is None or D_ait.shape[0] <= 800:
        pd.DataFrame(D_ait, index=sample_ids, columns=sample_ids).to_csv(
            OUT_DIR / f"dist_aitchison_matrix{suffix_full}.csv"
        )

    # PCoA
    X_ait, var_ait = cmdscale_pcoa(D_ait, k=10)
    X_br, var_br = cmdscale_pcoa(D_bray, k=10)
    pcoa_ait = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "PC1": X_ait[:, 0],
            "PC2": X_ait[:, 1],
            "pct_axis1_approx": var_ait[0] if len(var_ait) else np.nan,
            "pct_axis2_approx": var_ait[1] if len(var_ait) > 1 else np.nan,
        }
    )
    pcoa_ait.to_csv(OUT_DIR / f"pcoa_aitchison_coords{suffix_full}.tsv", sep="\t", index=False)
    pd.DataFrame(
        {"sample_id": sample_ids, "PC1": X_br[:, 0], "PC2": X_br[:, 1]}
    ).to_csv(OUT_DIR / f"pcoa_bray_coords{suffix_full}.tsv", sep="\t", index=False)

    # NMDS (very slow for large n; skipped by default in preview)
    stress_note = None
    nmds_ait = None
    if run_nmds:
        mds_ait = MDS(
            n_components=2,
            dissimilarity="precomputed",
            random_state=RNG_SEED,
            metric=False,
            n_init=2,
            max_iter=100,
            normalized_stress="auto",
        )
        nmds_ait = mds_ait.fit_transform(D_ait)
        stress_note = getattr(mds_ait, "stress_", None)
        with open(OUT_DIR / f"nmds_stress{suffix_full}.txt", "w", encoding="utf-8") as f:
            f.write(f"sklearn NMDS stress_ (Aitchison): {stress_note}\n")
        pd.DataFrame(
            {"sample_id": sample_ids, "NMDS1": nmds_ait[:, 0], "NMDS2": nmds_ait[:, 1]}
        ).to_csv(OUT_DIR / f"nmds_aitchison_scores{suffix_full}.tsv", sep="\t", index=False)
    else:
        with open(OUT_DIR / f"nmds_stress{suffix_full}.txt", "w", encoding="utf-8") as f:
            f.write("NMDS skipped (use full run or --nmds).\n")

    # PERMANOVA
    lines = []
    f1, p1, _ = permanova_stratified(D_ait, geo, study, geo_cats, n_perm, rng)
    lines.append(
        f"=== Model 1: {args.geo} | stratified permutations within Study (Aitchison) ===\n"
    )
    lines.append(f"pseudo-F = {f1:.6g}, p = {p1:.6g}, permutations = {n_perm}\n\n")

    f_study, p_study = permanova_unstratified(D_ait, study, study_cats, n_perm, rng)
    f_geo, p_geo = permanova_unstratified(D_ait, geo, geo_cats, n_perm, rng)
    lines.append("=== One-way PERMANOVA (unstratified, Aitchison) ===\n")
    lines.append(f"Study: pseudo-F = {f_study:.6g}, p = {p_study:.6g}\n")
    lines.append(f"{args.geo}: pseudo-F = {f_geo:.6g}, p = {p_geo:.6g}\n\n")
    lines.append(
        "Note: R adonis2 marginal order (Study+continent) differs from two separate one-way tests;\n"
        "use R script for exact Type-I partition if needed.\n"
    )
    with open(OUT_DIR / f"permanova_results_python{suffix_full}.txt", "w", encoding="utf-8") as f:
        f.writelines(lines)

    f_bd_a, p_bd_a = betadisper_perm(D_ait, geo, n_perm, rng)
    f_bd_b, p_bd_b = betadisper_perm(D_bray, geo, n_perm, rng)
    with open(OUT_DIR / f"permdisp_results_python{suffix_full}.txt", "w", encoding="utf-8") as f:
        f.write("=== betadisper-style: ANOVA on distance-to-centroid (PCoA r<=8), permuted F ===\n")
        f.write(f"group = {args.geo}\n")
        f.write(f"Aitchison: F = {f_bd_a:.6g}, p = {p_bd_a:.6g}\n")
        f.write(f"Bray-Curtis: F = {f_bd_b:.6g}, p = {p_bd_b:.6g}\n")

    # —— Dispersion: distance to centroid in PCoA space (betadisper-consistent r<=8) ——
    n_s = D_ait.shape[0]
    r_pc = min(8, n_s - 1)
    X_ait_r, _ = cmdscale_pcoa(D_ait, k=r_pc)
    X_br_r, _ = cmdscale_pcoa(D_bray, k=r_pc)
    disp_ait = _dispersion_to_centroid(X_ait_r, geo)
    disp_br = _dispersion_to_centroid(X_br_r, geo)
    pd.DataFrame(
        {
            "sample_id": sample_ids,
            "geo_group": geo,
            "dist_to_centroid_pcoa_aitchison": disp_ait,
            "dist_to_centroid_pcoa_bray": disp_br,
        }
    ).to_csv(OUT_DIR / f"dispersion_to_centroid{suffix_full}.tsv", sep="\t", index=False)

    # —— Plots: geographic coloring + centroids; within-group dispersion boxplots ——
    cats = list(geo_cats)
    n_cat = len(cats)
    cmap = qualitative_colors_n(n_cat)
    leg_fs = 7 if n_cat > 18 else (8 if n_cat > 12 else 10)

    def plot_pcoa_centroids(x1, x2, title, fname_base: str, var_pc: np.ndarray):
        fig_w = 10.8 if n_cat > 16 else (10.0 if n_cat > 12 else 9.0)
        if n_cat >= 7:
            fig_w = max(fig_w, 11.2 + min(n_cat, 30) * 0.06)
        fig, ax = plt.subplots(figsize=(fig_w, 6.8))
        pc1_pct = float(var_pc[0]) if len(var_pc) > 0 else 0.0
        pc2_pct = float(var_pc[1]) if len(var_pc) > 1 else 0.0
        ann_fs = max(7.0, min(10.5, 11.0 - n_cat / 5.5))
        for i, c in enumerate(cats):
            m = geo == c
            ax.scatter(
                x1[m],
                x2[m],
                s=20,
                alpha=0.55,
                label=str(c)[:40],
                color=cmap[i % len(cmap)],
            )
        cent_df = pd.DataFrame({"x": x1, "y": x2, "g": geo}).groupby("g", sort=False)[["x", "y"]].mean()
        cat_to_i = {c: j for j, c in enumerate(cats)}
        for gname, row in cent_df.iterrows():
            ci = cat_to_i.get(gname, 0)
            col = cmap[ci % len(cmap)]
            ax.scatter(
                row["x"],
                row["y"],
                s=240,
                marker="X",
                edgecolors="black",
                linewidths=0.9,
                facecolors=col,
                zorder=5,
            )

        arrow_kw = dict(
            arrowstyle="-",
            color="0.42",
            lw=0.75,
            shrinkA=3,
            shrinkB=1,
            connectionstyle="arc3,rad=0.0",
        )
        text_kw = dict(
            fontsize=ann_fs,
            fontweight="medium",
            color="#1a1a1a",
            bbox=dict(
                boxstyle="round,pad=0.22",
                facecolor="white",
                edgecolor="none",
                alpha=0.9,
            ),
            zorder=7,
        )

        # Centroid labels: >=7 groups listed vertically outside right edge + guide lines; fewer groups use short radial leaders to reduce overlap
        use_side_legend = n_cat >= 7
        items: list[tuple[object, float, float, int]] = []
        for c in cats:
            if c not in cent_df.index:
                continue
            row = cent_df.loc[c]
            items.append((c, float(row["x"]), float(row["y"]), cat_to_i[c]))

        ax.relim()
        ax.autoscale_view()
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
        xspan = xmax - xmin
        yspan = ymax - ymin
        span = max(xspan, yspan, 1e-9)

        if use_side_legend:
            trans_r = blended_transform_factory(ax.transAxes, ax.transData)
            # Sort centroids by y descending to match top-to-bottom text rows and reduce guide-line crossings
            items_sorted = sorted(items, key=lambda t: -t[2])
            n_lab = len(items_sorted)
            pad_y = 0.04 * yspan
            ys = np.linspace(ymax - pad_y, ymin + pad_y, n_lab) if n_lab > 1 else np.array([(ymax + ymin) / 2])
            for (gname, xc, yc, _), yt in zip(items_sorted, ys):
                gstr = str(gname).strip()
                lbl = gstr[:44] + ("\u2026" if len(gstr) > 44 else "")
                ax.annotate(
                    lbl,
                    xy=(xc, yc),
                    xytext=(1.01, yt),
                    textcoords=trans_r,
                    arrowprops=arrow_kw,
                    ha="left",
                    va="center",
                    annotation_clip=False,
                    **text_kw,
                )
        else:
            cx = float(cent_df["x"].mean())
            cy = float(cent_df["y"].mean())
            base_off = 0.055 * span
            n_it = len(items)
            jitter_deg = min(14.0, 55.0 / max(n_it, 1))
            for k, (gname, xc, yc, _) in enumerate(items):
                gstr = str(gname).strip()
                lbl = gstr[:44] + ("\u2026" if len(gstr) > 44 else "")
                dx, dy = xc - cx, yc - cy
                hn = float(np.hypot(dx, dy))
                if hn < 1e-12:
                    theta = 2 * np.pi * k / max(n_it, 1)
                    dx, dy = float(np.cos(theta)), float(np.sin(theta))
                else:
                    dx, dy = dx / hn, dy / hn
                deg = jitter_deg * (k - (n_it - 1) / 2.0)
                rad = np.deg2rad(deg)
                rdx = dx * np.cos(rad) - dy * np.sin(rad)
                rdy = dx * np.sin(rad) + dy * np.cos(rad)
                tx = xc + rdx * base_off * 2.8
                ty = yc + rdy * base_off * 2.8
                ax.annotate(
                    lbl,
                    xy=(xc, yc),
                    xytext=(tx, ty),
                    textcoords="data",
                    arrowprops=arrow_kw,
                    ha="center",
                    va="center",
                    annotation_clip=False,
                    **text_kw,
                )

        ax.set_xlabel(f"PCoA axis 1 ({pc1_pct:.1f}% variance explained)")
        ax.set_ylabel(f"PCoA axis 2 ({pc2_pct:.1f}% variance explained)")
        ax.set_title(title)
        if use_side_legend:
            ax.legend(
                bbox_to_anchor=(0.02, 0.98),
                loc="upper left",
                bbox_transform=ax.transAxes,
                fontsize=leg_fs,
                ncol=1,
                framealpha=0.92,
                fancybox=True,
            )
        else:
            ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=leg_fs, ncol=1)
        ax.grid(True, alpha=0.3)
        if use_side_legend:
            fig.tight_layout(rect=[0, 0, 0.74, 1])
        else:
            fig.tight_layout()
        for ext in (".png", ".pdf"):
            fig.savefig(OUT_DIR / f"{fname_base}{ext}", dpi=300 if ext == ".png" else None, bbox_inches="tight")
        plt.close(fig)

    def plot_dispersion_box(dist_vec: np.ndarray, title: str, fname_base: str):
        fig, ax = plt.subplots(figsize=(10, 5.5))
        data = [dist_vec[geo == c] for c in cats]
        bp = ax.boxplot(data, patch_artist=True, showfliers=False)
        ax.set_xticks(np.arange(1, len(cats) + 1))
        ax.set_xticklabels([str(c)[:25] for c in cats], rotation=45, ha="right", fontsize=10)
        for i, p in enumerate(bp["boxes"]):
            p.set_facecolor(cmap[i % len(cmap)])
            p.set_alpha(0.55)
        ax.set_ylabel("Distance to group centroid\n(PCoA space, r≤8)")
        ax.set_title(title)
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha="right", fontsize=10)
        ax.grid(True, axis="y", alpha=0.3)
        fig.tight_layout()
        for ext in (".png", ".pdf"):
            fig.savefig(OUT_DIR / f"{fname_base}{ext}", dpi=300 if ext == ".png" else None, bbox_inches="tight")
        plt.close(fig)

    pct1, pct2 = var_ait[0], var_ait[1] if len(var_ait) > 1 else 0.0
    prev_note = f"\n[preview: first {max_taxa} vOTUs]" if max_taxa is not None else ""
    geo_title = "Continent" if args.geo == "continent" else "Country"

    plot_pcoa_centroids(
        X_ait[:, 0],
        X_ait[:, 1],
        f"PCoA (Aitchison / CLR) by {geo_title}\nX = group centroid (labeled){prev_note}",
        f"fig_pcoa_aitchison_centroids_python{suffix_full}",
        var_ait,
    )
    plot_pcoa_centroids(
        X_br[:, 0],
        X_br[:, 1],
        f"PCoA (Bray-Curtis) by {geo_title}\nX = group centroid (labeled){prev_note}",
        f"fig_pcoa_bray_centroids_python{suffix_full}",
        var_br,
    )
    if not args.no_dual_pcoa:
        country_arr = meta["Country"].astype(str).str.strip().values
        cont_arr = meta["continent"].astype(str).values
        k_plot = int(args.country_plot_top)
        co_lbl = "all countries" if k_plot <= 0 else f"top-{k_plot}+Other"
        plot_pcoa_country_color_continent_shape(
            X_ait[:, 0],
            X_ait[:, 1],
            country_arr,
            cont_arr,
            f"PCoA (Aitchison) | color=Country ({co_lbl}), marker=Continent\n"
            f"~{pct1:.1f}% / ~{pct2:.1f}% var. (approx.){prev_note}",
            OUT_DIR / f"fig_pcoa_aitchison_countryColor_continentShape{suffix_full}",
            k_plot,
        )
        plot_pcoa_country_color_continent_shape(
            X_br[:, 0],
            X_br[:, 1],
            country_arr,
            cont_arr,
            f"PCoA (Bray-Curtis) | color=Country ({co_lbl}), marker=Continent{prev_note}",
            OUT_DIR / f"fig_pcoa_bray_countryColor_continentShape{suffix_full}",
            k_plot,
        )
    plot_dispersion_box(
        disp_ait,
        f"Within-group dispersion ({geo_title}): distance to centroid\n"
        f"Aitchison PCoA (PERMDISP-style; F={f_bd_a:.3f}, p={p_bd_a:.4f}){prev_note}",
        f"fig_dispersion_aitchison_boxplot_python{suffix_full}",
    )
    plot_dispersion_box(
        disp_br,
        f"Within-group dispersion ({geo_title}): distance to centroid\n"
        f"Bray PCoA (PERMDISP-style; F={f_bd_b:.3f}, p={p_bd_b:.4f}){prev_note}",
        f"fig_dispersion_bray_boxplot_python{suffix_full}",
    )

    # 2x2 combined panel (for submission)
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    for ax, (x1, x2, tit) in zip(
        axes[0],
        [
            (X_ait[:, 0], X_ait[:, 1], f"Aitchison PCoA ({geo_title})"),
            (X_br[:, 0], X_br[:, 1], f"Bray PCoA ({geo_title})"),
        ],
    ):
        for i, c in enumerate(cats):
            m = geo == c
            ax.scatter(x1[m], x2[m], s=14, alpha=0.5, color=cmap[i % len(cmap)], label=str(c)[:20])
        cent = pd.DataFrame({"x": x1, "y": x2, "g": geo}).groupby("g")[["x", "y"]].mean()
        cti = {c: j for j, c in enumerate(cats)}
        for gname, row in cent.iterrows():
            idx = cti.get(gname, 0)
            ax.scatter(
                row["x"],
                row["y"],
                s=120,
                marker="X",
                edgecolors="k",
                linewidths=0.6,
                facecolors=cmap[idx % len(cmap)],
                zorder=5,
            )
        ax.set_title(tit)
        ax.set_xlabel("Axis 1")
        ax.set_ylabel("Axis 2")
        ax.grid(True, alpha=0.25)
    axes[0, 0].legend(fontsize=7, ncol=2, loc="upper left", bbox_to_anchor=(0, 1.25))
    for ax, dist_v, tit in zip(
        axes[1],
        [disp_ait, disp_br],
        ["Dispersion (Aitchison PCoA)", "Dispersion (Bray PCoA)"],
    ):
        data = [dist_v[geo == c] for c in cats]
        bp = ax.boxplot(data, patch_artist=True, showfliers=False)
        ax.set_xticks(np.arange(1, len(cats) + 1))
        ax.set_xticklabels([str(c)[:12] for c in cats], rotation=60, ha="right", fontsize=9)
        for i, p in enumerate(bp["boxes"]):
            p.set_facecolor(cmap[i % len(cmap)])
            p.set_alpha(0.5)
        ax.set_ylabel("Dist. to centroid")
        ax.set_title(tit)
        ax.grid(True, axis="y", alpha=0.3)
    fig.suptitle(
        f"β diversity: centroid position & within-group dispersion ({geo_title})",
        fontsize=13,
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(OUT_DIR / f"fig_beta_diversity_panel_python{suffix_full}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / f"fig_beta_diversity_panel_python{suffix_full}.pdf", bbox_inches="tight")
    plt.close(fig)

    if nmds_ait is not None:
        fig, ax = plt.subplots(figsize=(8.5, 6))
        for i, c in enumerate(cats):
            m = geo == c
            ax.scatter(
                nmds_ait[m, 0],
                nmds_ait[m, 1],
                s=18,
                alpha=0.6,
                label=str(c)[:35],
                color=cmap[i % len(cmap)],
            )
        ax.set_title(f"NMDS (Aitchison); stress={stress_note}{prev_note}")
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=leg_fs)
        fig.tight_layout()
        fig.savefig(
            OUT_DIR / f"fig_nmds_aitchison_python{suffix_full}.png",
            dpi=300,
            bbox_inches="tight",
        )
        fig.savefig(OUT_DIR / f"fig_nmds_aitchison_python{suffix_full}.pdf", bbox_inches="tight")
        plt.close(fig)

    with open(OUT_DIR / f"interpretation_summary_python{suffix_full}.txt", "w", encoding="utf-8") as f:
        f.write(
            f"Geo grouping: {args.geo}. "
            "Aitchison = Euclidean on CLR (+pseudocount). "
            "Stratified PERMANOVA shuffles geo labels within each Study.\n"
            "Dispersion = distance to group centroid in first r≤8 PCoA axes (betadisper-style).\n"
        )

    print("Done. Output directory:", OUT_DIR, flush=True)


if __name__ == "__main__":
    main()
