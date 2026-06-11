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
Distance–decay (geographic vs community dissimilarity) for selected bioprojects.

- PRJEB11755: multinational → full study + Denmark-only + France-only + cross-country pair plot
  (China arm has only one Region in meta → no within-China regional decay.)
- CNP0005498, PRJNA684454, PRJEB62878: single-country multi-region → within-study decay

Geographic distance: great-circle km from region_latlon_lookup.tsv (approximate centroids).
Community: Aitchison (CLR + Euclidean) and Bray–Curtis on relative abundances.
Mantel: Spearman correlation of upper-triangle distances, label permutation of samples.

Usage:
  python distance_decay_spatial_subsets.py
  python distance_decay_spatial_subsets.py --max-taxa 10000 --perm 999

Env: PHAGE_PROJECT_ROOT (Phage folder)
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import pairwise_distances

# -----------------------------------------------------------------------------
# paths
# -----------------------------------------------------------------------------
_THIS = Path(__file__).resolve()
_ROOT = Path(os.environ.get("PHAGE_PROJECT_ROOT", str(_THIS.parent.parent))).expanduser().resolve()
MATRIX = _ROOT / "FINAL_all_projects_Count_matrix.tsv"
META = _ROOT / "pvca_sample_meta_full.tsv"
LOOKUP = _THIS.parent / "region_latlon_lookup.tsv"
OUT = _THIS.parent / "distance_decay_spatial"
PSEUDOCOUNT = 1.0
RNG_SEED = 42


def _strip_count_col(name: str) -> str:
    return name.replace("_Count", "") if name.endswith("_Count") else name


def load_latlon(path: Path) -> dict[str, tuple[float, float]]:
    rows = []
    with path.open(encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3 or parts[0] == "Region":
                continue
            reg, lat_s, lon_s = parts[0], parts[1], parts[2]
            rows.append((reg, float(lat_s), float(lon_s)))
    return {r: (la, lo) for r, la, lo in rows}


def _matrix_col_indices(matrix_path: Path, sample_ids: list[str]) -> tuple[list[str], list[int]]:
    with matrix_path.open(encoding="utf-8", errors="replace") as f:
        header = f.readline().strip().split("\t")
    if not header or header[0] != "Contig":
        raise ValueError("Matrix must have Contig as first column")
    id_to_idx = {_strip_count_col(h): i for i, h in enumerate(header[1:], start=1)}
    miss = [s for s in sample_ids if s not in id_to_idx]
    if miss:
        raise ValueError(f"Samples not in matrix header: {miss[:5]}… ({len(miss)} total)")
    idxs = [id_to_idx[s] for s in sample_ids if s in id_to_idx]
    return header, idxs


def read_matrix_for_samples(
    matrix_path: Path, sample_ids: list[str], max_taxa: int
) -> tuple[np.ndarray, list[str], list[str]]:
    """
    Two-pass line-wise read: (1) total counts per contig across selected samples,
    (2) keep top max_taxa contigs and build samples x taxa matrix.
    """
    header, col_idx = _matrix_col_indices(matrix_path, sample_ids)
    n_samp = len(col_idx)
    sums: dict[str, float] = {}
    with matrix_path.open(encoding="utf-8", errors="replace") as f:
        f.readline()
        for line in f:
            if not line.strip():
                continue
            cells = line.rstrip("\n").split("\t")
            if len(cells) < max(col_idx) + 1:
                continue
            cid = cells[0]
            tot = 0.0
            for j in col_idx:
                try:
                    tot += float(cells[j])
                except ValueError:
                    pass
            sums[cid] = tot
    if not sums:
        raise ValueError("No matrix rows read")
    ordered = sorted(sums.keys(), key=lambda k: sums[k], reverse=True)[:max_taxa]
    want = set(ordered)
    rows_list: list[list[float]] = []
    taxa: list[str] = []
    with matrix_path.open(encoding="utf-8", errors="replace") as f:
        f.readline()
        for line in f:
            if not line.strip():
                continue
            cells = line.rstrip("\n").split("\t")
            cid = cells[0]
            if cid not in want:
                continue
            row = []
            ok = True
            for j in col_idx:
                try:
                    row.append(float(cells[j]))
                except (ValueError, IndexError):
                    row.append(0.0)
            rows_list.append(row)
            taxa.append(cid)
    # enforce same order as ordered (top abundance)
    pos = {t: i for i, t in enumerate(taxa)}
    order_idx = [pos[t] for t in ordered if t in pos]
    mat = np.array([rows_list[i] for i in order_idx], dtype=np.float64).T
    taxa_ordered = [taxa[i] for i in order_idx]
    matched_ids = [_strip_count_col(header[j]) for j in col_idx]
    return mat, matched_ids, taxa_ordered


def clr_rows(Y: np.ndarray, pseudocount: float = PSEUDOCOUNT) -> np.ndarray:
    z = np.log(Y + pseudocount)
    return z - z.mean(axis=1, keepdims=True)


def haversine_km(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    """Pairwise great-circle km, shape (n,n)."""
    la = np.radians(lat)
    lo = np.radians(lon)
    n = len(lat)
    out = np.zeros((n, n), dtype=np.float64)
    r = 6371.0
    for i in range(n - 1):
        dlat = la[i + 1 :] - la[i]
        dlon = lo[i + 1 :] - lo[i]
        a = np.sin(dlat / 2) ** 2 + np.cos(la[i]) * np.cos(la[i + 1 :]) * np.sin(dlon / 2) ** 2
        a = np.minimum(1.0, a)
        d = 2 * r * np.arcsin(np.sqrt(a))
        out[i, i + 1 :] = d
        out[i + 1 :, i] = d
    return out


def mantel_spearman(d1: np.ndarray, d2: np.ndarray, perm: int, rng: np.random.Generator) -> tuple[float, float]:
    n = d1.shape[0]
    tri = np.triu_indices(n, k=1)
    x, y = d1[tri], d2[tri]
    r_obs = float(spearmanr(x, y, nan_policy="omit").correlation)
    if np.isnan(r_obs):
        return float("nan"), float("nan")
    ge = 0
    for _ in range(perm):
        p = rng.permutation(n)
        dp = d1[np.ix_(p, p)]
        rp = spearmanr(dp[tri], y, nan_policy="omit").correlation
        if np.isnan(rp):
            continue
        if abs(rp) >= abs(r_obs):
            ge += 1
    pval = (1 + ge) / (1 + perm)
    return r_obs, pval


def partial_mantel_spearman(
    d_comm: np.ndarray, d_geo: np.ndarray, d_ctrl: np.ndarray, perm: int, rng: np.random.Generator
) -> tuple[float, float]:
    """Residual Spearman correlation after rank-based partialing (approximate; for depth control)."""
    n = d_comm.shape[0]
    tri = np.triu_indices(n, k=1)
    from scipy.stats import rankdata

    def rvec(d):
        v = d[tri]
        return rankdata(v)

    rc, rg, rz = rvec(d_comm), rvec(d_geo), rvec(d_ctrl)
    # partial corr via linear regression on ranks (approximate partial Mantel)
    X = np.column_stack([np.ones(len(rz)), rz])
    beta_c, _, _, _ = np.linalg.lstsq(X, rc, rcond=None)
    beta_g, _, _, _ = np.linalg.lstsq(X, rg, rcond=None)
    rc_r = rc - X @ beta_c
    rg_r = rg - X @ beta_g
    r_obs = float(spearmanr(rc_r, rg_r, nan_policy="omit").correlation)
    if np.isnan(r_obs):
        return float("nan"), float("nan")
    ge = 0
    for _ in range(perm):
        p = rng.permutation(n)
        dp = d_comm[np.ix_(p, p)]
        rc2 = rankdata(dp[tri])
        beta_c2, _, _, _ = np.linalg.lstsq(X, rc2, rcond=None)
        rc2r = rc2 - X @ beta_c2
        rp = spearmanr(rc2r, rg_r, nan_policy="omit").correlation
        if np.isnan(rp):
            continue
        if abs(rp) >= abs(r_obs):
            ge += 1
    return r_obs, (1 + ge) / (1 + perm)


def run_one_subset(
    label: str,
    meta_sub: pd.DataFrame,
    latlon: dict[str, tuple[float, float]],
    max_taxa: int,
    perm: int,
    rng: np.random.Generator,
) -> dict:
    sids = meta_sub["sample_id"].astype(str).tolist()
    Y, matched_ids, taxa = read_matrix_for_samples(MATRIX, sids, max_taxa)
    meta_sub = meta_sub.set_index("sample_id").loc[matched_ids].reset_index()

    regions = meta_sub["Region"].astype(str).tolist()
    missing = [r for r in set(regions) if r not in latlon]
    if missing:
        raise ValueError(f"{label}: missing lat/lon for regions: {missing[:8]}")

    lat = np.array([latlon[r][0] for r in regions], dtype=np.float64)
    lon = np.array([latlon[r][1] for r in regions], dtype=np.float64)
    d_geo = haversine_km(lat, lon)
    depths = Y.sum(axis=1)
    logd = np.log(depths + 1.0)
    d_depth = np.abs(np.subtract.outer(logd, logd))

    clr = clr_rows(Y)
    d_ait = squareform(pdist(clr, metric="euclidean"))
    rel = Y / np.maximum(Y.sum(axis=1, keepdims=True), 1e-12)
    d_bray = pairwise_distances(rel, metric="braycurtis", n_jobs=1)

    tri = np.triu_indices(len(matched_ids), k=1)
    g = d_geo[tri]
    ba = d_bray[tri]
    aa = d_ait[tri]

    r_m_ait, p_m_ait = mantel_spearman(d_ait, d_geo, perm, rng)
    r_m_br, p_m_br = mantel_spearman(d_bray, d_geo, perm, rng)
    r_pm_ait, p_pm_ait = partial_mantel_spearman(d_ait, d_geo, d_depth, perm, rng)
    r_pm_br, p_pm_br = partial_mantel_spearman(d_bray, d_geo, d_depth, perm, rng)

    return {
        "label": label,
        "n_samples": len(matched_ids),
        "n_taxa": Y.shape[1],
        "mantel_ait_r": r_m_ait,
        "mantel_ait_p": p_m_ait,
        "mantel_bray_r": r_m_br,
        "mantel_bray_p": p_m_br,
        "partial_mantel_ait_r": r_pm_ait,
        "partial_mantel_ait_p": p_pm_ait,
        "partial_mantel_bray_r": r_pm_br,
        "partial_mantel_bray_p": p_pm_br,
        "plot_geo": g,
        "plot_bray": ba,
        "plot_ait": aa,
        "d_geo_full": d_geo,
        "d_bray_full": d_bray,
        "d_ait_full": d_ait,
    }


def plot_hex(geo, comm, title: str, out_stem: Path, xlab: str = "Geographic distance (km)", ylab: str = ""):
    """Save PNG + PDF (same stem). out_stem = path without suffix, e.g. OUT / decay_foo_bray"""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42})

    max_pts = 80000
    if len(geo) > max_pts:
        idx = np.random.default_rng(0).choice(len(geo), size=max_pts, replace=False)
        geo, comm = geo[idx], comm[idx]
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    hb = ax.hexbin(geo, comm, gridsize=45, mincnt=1, cmap="viridis")
    plt.colorbar(hb, ax=ax, label="count")
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(out_stem.with_suffix(".png"), dpi=200, bbox_inches="tight")
    fig.savefig(out_stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-taxa", type=int, default=10000, help="Keep top N vOTUs by total counts")
    ap.add_argument("--perm", type=int, default=999, help="Mantel permutations")
    args = ap.parse_args()
    rng = np.random.default_rng(RNG_SEED)

    OUT.mkdir(parents=True, exist_ok=True)
    latlon = load_latlon(LOOKUP)
    meta = pd.read_csv(META, sep="\t", dtype=str)
    for c in meta.columns:
        meta[c] = meta[c].fillna("").str.strip()

    rows_out: list[dict] = []

    configs: list[tuple[str, pd.DataFrame, str]] = []

    # PRJEB11755
    m = meta[meta["Project"] == "PRJEB11755"].copy()
    configs.append(("PRJEB11755_all", m, ""))
    md = m[m["Country"] == "Denmark"]
    configs.append(("PRJEB11755_Denmark", md, ""))
    mf = m[m["Country"] == "France"]
    configs.append(("PRJEB11755_France", mf, ""))

    for proj in ["CNP0005498", "PRJNA684454", "PRJEB62878"]:
        configs.append((f"{proj}_within_country", meta[meta["Project"] == proj].copy(), ""))

    def save_plots(res: dict, label: str) -> None:
        suf = label.replace(" ", "_")
        plot_hex(
            res["plot_geo"],
            res["plot_bray"],
            f"{label}\nBray–Curtis vs geographic km\nMantel r={res['mantel_bray_r']:.3f}, p={res['mantel_bray_p']:.4f}",
            OUT / f"decay_{suf}_bray",
            ylab="Bray–Curtis",
        )
        plot_hex(
            res["plot_geo"],
            res["plot_ait"],
            f"{label}\nAitchison vs geographic km\nMantel r={res['mantel_ait_r']:.3f}, p={res['mantel_ait_p']:.4f}",
            OUT / f"decay_{suf}_aitchison",
            ylab="Aitchison (Euclidean on CLR)",
        )

    scalar_keys = [
        "label",
        "n_samples",
        "n_taxa",
        "mantel_ait_r",
        "mantel_ait_p",
        "mantel_bray_r",
        "mantel_bray_p",
        "partial_mantel_ait_r",
        "partial_mantel_ait_p",
        "partial_mantel_bray_r",
        "partial_mantel_bray_p",
    ]

    for label, msub, _ in configs:
        if len(msub) < 6:
            print(f"SKIP {label}: n={len(msub)}", flush=True)
            continue
        try:
            res = run_one_subset(label, msub, latlon, args.max_taxa, args.perm, rng)
        except Exception as e:
            print(f"FAIL {label}: {e}", flush=True)
            continue

        row = {k: res[k] for k in scalar_keys}
        row["test_type"] = "mantel_full_matrix"
        rows_out.append(row)
        save_plots(res, label)
        print(f"OK {label} n={res['n_samples']} taxa={res['n_taxa']}", flush=True)

        # Cross-border pairs only (descriptive Spearman; full-sample Mantel is in PRJEB11755_all)
        if label == "PRJEB11755_all":
            countries = msub["Country"].astype(str).values
            tri = np.triu_indices(res["n_samples"], k=1)
            mask = countries[tri[0]] != countries[tri[1]]
            g = res["d_geo_full"][tri][mask]
            ba = res["d_bray_full"][tri][mask]
            aa = res["d_ait_full"][tri][mask]
            r_b, p_b = spearmanr(g, ba, nan_policy="omit")
            r_a, p_a = spearmanr(g, aa, nan_policy="omit")
            rows_out.append(
                {
                    "label": "PRJEB11755_cross_country_pairs_only",
                    "test_type": "spearman_pair_subset",
                    "n_samples": res["n_samples"],
                    "n_taxa": res["n_taxa"],
                    "n_pairs_cross_country": int(mask.sum()),
                    "mantel_ait_r": float(r_a) if not np.isnan(r_a) else np.nan,
                    "mantel_ait_p": float(p_a) if not np.isnan(p_a) else np.nan,
                    "mantel_bray_r": float(r_b) if not np.isnan(r_b) else np.nan,
                    "mantel_bray_p": float(p_b) if not np.isnan(p_b) else np.nan,
                    "partial_mantel_ait_r": np.nan,
                    "partial_mantel_ait_p": np.nan,
                    "partial_mantel_bray_r": np.nan,
                    "partial_mantel_bray_p": np.nan,
                    "note": "Spearman(geo, beta) on i<j with Country_i!=Country_j; not Mantel",
                }
            )
            plot_hex(
                g,
                ba,
                "PRJEB11755 cross-country pairs\nBray–Curtis vs km\n"
                f"Spearman r={r_b:.3f}, p={p_b:.4f} (pair subset)",
                OUT / "decay_PRJEB11755_crosscountry_bray",
                ylab="Bray–Curtis",
            )
            plot_hex(
                g,
                aa,
                "PRJEB11755 cross-country pairs\nAitchison vs km\n"
                f"Spearman r={r_a:.3f}, p={p_a:.4f} (pair subset)",
                OUT / "decay_PRJEB11755_crosscountry_aitchison",
                ylab="Aitchison (Euclidean on CLR)",
            )

    cols_order = [
        "label",
        "test_type",
        "n_samples",
        "n_taxa",
        "n_pairs_cross_country",
        "mantel_ait_r",
        "mantel_ait_p",
        "mantel_bray_r",
        "mantel_bray_p",
        "partial_mantel_ait_r",
        "partial_mantel_ait_p",
        "partial_mantel_bray_r",
        "partial_mantel_bray_p",
        "note",
    ]
    df_sum = pd.DataFrame(rows_out)
    for c in cols_order:
        if c not in df_sum.columns:
            df_sum[c] = np.nan
    df_sum[cols_order].to_csv(OUT / "distance_decay_mantel_summary.tsv", sep="\t", index=False)

    note = OUT / "README_distance_decay.txt"
    note.write_text(
        "Geographic coordinates are approximate region/city centroids (see region_latlon_lookup.tsv).\n"
        "Partial Mantel here = rank-residual approach controlling log(sequencing depth); "
        "for strict vegan::mantel.partial replicate in R if needed.\n"
        "PRJEB11755 China: only one Region in metadata → no within-China regional decay for this study.\n",
        encoding="utf-8",
    )
    print(f"Wrote {OUT / 'distance_decay_mantel_summary.tsv'}", flush=True)


if __name__ == "__main__":
    main()
