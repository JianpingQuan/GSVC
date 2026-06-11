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
Sample-based vOTU accumulation curves by country (sample-based accumulation / rarefaction-style)

This script, output figures, and TSVs live in 03_Biogeographical/ for easy comparison with results.
Matrix and biogeography_descriptive_analysis.py are in the parent Phage directory; paths resolved via __file__.

Method: accumulate samples in random order per country, merge vOTUs with count>0 into a set and record union size; repeat many rounds for mean/quantiles.
Note: accumulation by sample count, not sequencing read-depth rarefaction.

Outputs (this directory):
  - biogeo_rarefaction_by_country_curve.tsv
  - biogeo_rarefaction_by_country_summary.tsv
  - biogeo_fig_rarefaction_by_country.pdf / .png

With --max-rows N, filenames append _preview{N}.

Run (from this directory):
  python biogeography_rarefaction_by_country.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_SCRIPT_DIR = Path(__file__).resolve().parent
OUT_DIR = _SCRIPT_DIR
BASE = _SCRIPT_DIR.parent
sys.path.insert(0, str(BASE))
from biogeography_descriptive_analysis import (  # noqa: E402
    BIOGEO_COUNTRY_PALETTE,
    MATRIX,
    load_meta,
    safe_unlink,
)

CHUNKSIZE = 250
N_PERM = 80
RNG_SEED = 42


def read_matrix_header(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
    if not header or header[0] != "Contig":
        raise ValueError("First column should be Contig")
    return header


def save_figure_dual(fig, stem: Path) -> None:
    import matplotlib.pyplot as plt

    stem.parent.mkdir(parents=True, exist_ok=True)
    p_pdf = stem.with_suffix(".pdf")
    p_png = stem.with_suffix(".png")
    safe_unlink(p_pdf)
    safe_unlink(p_png)
    try:
        fig.savefig(p_pdf, bbox_inches="tight", format="pdf", facecolor="white", edgecolor="none")
        fig.savefig(p_png, bbox_inches="tight", format="png", facecolor="white", edgecolor="none", dpi=300)
    except OSError:
        alt = stem.parent / "_write_ok"
        alt.mkdir(parents=True, exist_ok=True)
        p_pdf = alt / (stem.name + ".pdf")
        p_png = alt / (stem.name + ".png")
        fig.savefig(p_pdf, bbox_inches="tight", format="pdf", facecolor="white", edgecolor="none")
        fig.savefig(p_png, bbox_inches="tight", format="png", facecolor="white", edgecolor="none", dpi=300)
    plt.close(fig)
    print("Writed:", p_pdf)
    print("Writed:", p_png)


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    df.to_csv(tmp, sep="\t", index=False)
    safe_unlink(path)
    try:
        tmp.replace(path)
        return
    except OSError:
        pass
    try:
        safe_unlink(path)
        tmp.replace(path)
        return
    except OSError:
        pass
    alt_dir = path.parent / "_write_ok"
    alt_dir.mkdir(parents=True, exist_ok=True)
    alt = alt_dir / path.name
    try:
        tmp.replace(alt)
    except OSError:
        df.to_csv(alt, sep="\t", index=False)
        try:
            tmp.unlink(missing_ok=True)
        except OSError:
            pass


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot sample-based cumulative vOTU curves by country (all countries in the metadata)"
    )
    parser.add_argument(
        "--max-rows",
        type=int,
        default=None,
        metavar="N",
        help="Read only the first N rows of the vOTU matrix (quick preview); output filename includes the _preview{N} suffix.",
    )
    args = parser.parse_args()
    max_rows = args.max_rows
    out_suffix = f"_preview{max_rows}" if max_rows is not None else ""

    from biogeography_descriptive_analysis import apply_nature_style

    import matplotlib.pyplot as plt

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    apply_nature_style()

    meta = load_meta()
    header = read_matrix_header(MATRIX)
    sample_ids = [
        c[: -len("_Count")] if c.endswith("_Count") else c for c in header[1:]
    ]
    n_samples = len(sample_ids)
    sid_to_idx = {s: i for i, s in enumerate(sample_ids)}

    meta = meta.drop_duplicates("sample_id")
    country_by_idx = ["NA"] * n_samples
    for row in meta.itertuples(index=False):
        sid = str(row.sample_id)
        cty = str(row.Country) if getattr(row, "Country", None) is not None else "NA"
        j = sid_to_idx.get(sid)
        if j is not None:
            country_by_idx[j] = cty

    all_meta_countries = sorted({c for c in country_by_idx if c and c != "NA"})
    country_to_indices: dict[str, list[int]] = {c: [] for c in all_meta_countries}
    for j, c in enumerate(country_by_idx):
        if c != "NA":
            country_to_indices[c].append(j)

    count_cols = header[1:]
    dtype_counts = {c: np.int32 for c in count_cols}

    if max_rows is not None:
        print(f"Constructing per-sample vOTU sets (streaming, previewing only the first {max_rows} rows)…", flush=True)
    else:
        print("Constructing per-sample vOTU sets (streaming, full matrix)…", flush=True)
    sample_sets: list[set[int]] = [set() for _ in range(n_samples)]
    row_id = 0
    chunk_i = 0
    for chunk in pd.read_csv(
        MATRIX,
        sep="\t",
        chunksize=CHUNKSIZE,
        dtype=dtype_counts,
        converters={"Contig": str},
        low_memory=False,
    ):
        data = chunk.drop(columns=["Contig"]).to_numpy(dtype=np.int32, copy=False)
        n_take = data.shape[0]
        if max_rows is not None:
            left = max_rows - row_id
            if left <= 0:
                break
            n_take = min(n_take, left)
        for r in range(n_take):
            nz = np.flatnonzero(data[r] > 0)
            rid = row_id + r
            for s in nz:
                sample_sets[int(s)].add(rid)
        row_id += n_take
        chunk_i += 1
        if chunk_i % 100 == 0:
            print(f"  … {row_id} 行", flush=True)
        if max_rows is not None and row_id >= max_rows:
            break

    countries = all_meta_countries
    rng = np.random.default_rng(RNG_SEED)
    curve_rows = []
    summary_rows = []

    print("Calculate rarefaction/accumulation curves for each country …", flush=True)
    for cty in countries:
        idx = country_to_indices.get(cty, [])
        k = len(idx)
        if k == 0:
            continue
        arr = np.array(idx, dtype=np.int32)
        total_union: set[int] = set()
        for j in idx:
            total_union |= sample_sets[j]
        total_richness = len(total_union)

        if k == 1:
            curves = np.array([[total_richness]], dtype=np.float64)
        else:
            n_perm = max(20, min(N_PERM, max(1, 80000 // k)))
            curves = np.zeros((n_perm, k), dtype=np.float64)
            for p in range(n_perm):
                perm = rng.permutation(arr)
                u: set[int] = set()
                for step, j in enumerate(perm):
                    u |= sample_sets[int(j)]
                    curves[p, step] = len(u)

        mean_c = curves.mean(axis=0)
        std_c = curves.std(axis=0)
        p_low = np.percentile(curves, 2.5, axis=0)
        p_high = np.percentile(curves, 97.5, axis=0)
        n_drawn = np.arange(1, k + 1)

        if k >= 10:
            w = max(1, k // 10)
            slope_end = (mean_c[-1] - mean_c[-w - 1]) / w if w < k else np.nan
        else:
            slope_end = (mean_c[-1] - mean_c[0]) / (k - 1) if k > 1 else 0.0

        summary_rows.append(
            {
                "Country": cty,
                "n_samples": k,
                "total_vOTUs_union": total_richness,
                "mean_richness_at_max_samples": float(mean_c[-1]),
                "approx_new_vOTU_per_sample_last_decile": float(slope_end),
            }
        )

        for a, b, lo, hi, s in zip(n_drawn, mean_c, p_low, p_high, std_c):
            curve_rows.append(
                {
                    "Country": cty,
                    "n_samples_accumulated": int(a),
                    "mean_observed_vOTUs": float(b),
                    "sd": float(s),
                    "p025": float(lo),
                    "p975": float(hi),
                }
            )

    curve_df = pd.DataFrame(curve_rows)
    summ_df = pd.DataFrame(summary_rows).sort_values("n_samples", ascending=False)

    curve_tsv = OUT_DIR / f"biogeo_rarefaction_by_country_curve{out_suffix}.tsv"
    summ_tsv = OUT_DIR / f"biogeo_rarefaction_by_country_summary{out_suffix}.tsv"
    write_tsv(curve_df, curve_tsv)
    write_tsv(summ_df, summ_tsv)
    print(f"已写: {curve_tsv.name}", flush=True)
    print(f"已写: {summ_tsv.name}", flush=True)

    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    palette = BIOGEO_COUNTRY_PALETTE
    plot_order = [c for c in countries if c in set(summ_df["Country"])]
    for i, cty in enumerate(plot_order):
        sub = curve_df[curve_df["Country"] == cty]
        col = palette[i % len(palette)]
        ns = int(summ_df.loc[summ_df["Country"] == cty, "n_samples"].iloc[0])
        ax.plot(
            sub["n_samples_accumulated"],
            sub["mean_observed_vOTUs"],
            color=col,
            linewidth=1.2,
            label=f"{cty} (n={ns})",
            alpha=0.9,
        )
    ax.set_xlabel("Number of samples (random accumulation order)")
    ax.set_ylabel("Observed vOTUs (cumulative union)")
    title = (
        "Sample-based vOTU accumulation curves by country\n"
        "(mean of random sample orders; not read-depth rarefaction)"
    )
    if max_rows is not None:
        title += f"\n[preview: first {max_rows} vOTU rows only]"
    ax.set_title(title)
    ax.grid(True, axis="y", linestyle="-", alpha=0.35)
    leg_fs = 8 if len(summ_df) <= 5 else 6
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=leg_fs, ncol=1)
    fig.tight_layout()
    fig_stem = OUT_DIR / f"biogeo_fig_rarefaction_by_country{out_suffix}"
    save_figure_dual(fig, fig_stem)

    print("Finished", flush=True)


if __name__ == "__main__":
    main()
