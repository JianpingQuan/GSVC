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
Descriptive biogeography analysis (streaming read of Count matrix)
- Sampling intensity, alpha diversity (observed, Shannon)
- Phylum/class/order/family/genus: stacked bar charts of country-mean sample composition
- Figures: Nature style, PDF + PNG output (300 dpi)

Output directory: 03_Biogeographical/
"""
from __future__ import annotations

import os
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd


def _project_root() -> Path:
    """Project root (contains matrix and meta). Prefer PHAGE_PROJECT_ROOT env var, else this file's directory."""
    env = os.environ.get("PHAGE_PROJECT_ROOT", "").strip()
    if env:
        return Path(env).expanduser().resolve()
    return Path(__file__).resolve().parent


BASE = _project_root()
MATRIX = BASE / "FINAL_all_projects_Count_matrix.tsv"
VIRUS_SUMMARY = BASE / "high_medium_quality_viral_rep_seq_virus_summary.tsv"
META_XLSX = BASE / "analysis_id_matched2.xlsx"
META_TSV = BASE / "pvca_sample_meta_full.tsv"
OUT_DIR = BASE / "03_Biogeographical"

CHUNKSIZE = 250

RANK_KEYS = ("phylum", "class", "order", "family", "genus")
RANK_LABELS_EN = {
    "phylum": "Phylum",
    "class": "Class",
    "order": "Order",
    "family": "Family",
    "genus": "Genus",
}
# Stacked chart: show Top N per level, merge rest into Other (increase Other at genus if annotations are sparse)
TOP_N_BY_RANK = {"phylum": 12, "class": 14, "order": 12, "family": 16, "genus": 16}

# Tableau / Nature-friendly low-saturation palette (cyclable)
NATURE_PALETTE = [
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
    "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
    "#499894", "#E377C2", "#8C564B", "#BCBD22", "#17BECF",
    "#9467BD", "#7F7F7F", "#D62728", "#1F77B4", "#2CA02C",
]

# Shared country-level curves / boxplots (consistent with biogeography_rarefaction_by_country.py; colors by country name alphabetical order)
BIOGEO_COUNTRY_PALETTE = [
    "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
    "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
    "#499894", "#E377C2", "#8C564B", "#BCBD22", "#17BECF",
    "#9467BD", "#D62728", "#1F77B4", "#2CA02C", "#FF7F0E",
]


def country_color_map(country_list: list) -> dict[str, str]:
    """Country -> hex color. Index by country name alphabetical order, consistent with rarefaction plot."""
    uniq = sorted(
        {
            str(c)
            for c in country_list
            if c is not None and str(c) != "NA" and str(c).strip() != ""
        }
    )
    pal = BIOGEO_COUNTRY_PALETTE
    return {c: pal[i % len(pal)] for i, c in enumerate(uniq)}


def lighten_hex_color(hex_color: str, white_mix: float = 0.52) -> str:
    """Blend with white for boxplot fill (outline keeps original color)."""
    import matplotlib.colors as mcolors

    r, g, b = mcolors.to_rgb(hex_color)
    w = white_mix
    return mcolors.to_hex((w + (1 - w) * r, w + (1 - w) * g, w + (1 - w) * b))


def _default_ranks() -> dict[str, str]:
    return {k: "Unclassified" for k in RANK_KEYS}


def parse_all_ranks(tax: str) -> dict[str, str]:
    """Parse phylum/class/order/family from semicolon-separated ICTV-style taxonomy; genus: terminal *virus but not *viridae."""
    out = _default_ranks()
    if pd.isna(tax) or not str(tax).strip():
        return out
    parts = [x.strip() for x in str(tax).split(";") if x.strip()]
    if not parts:
        return out
    for p in parts:
        if re.search(r"viricota$", p, re.I):
            out["phylum"] = p
        if re.search(r"viricetes$", p, re.I):
            out["class"] = p
        if re.search(r"virales$", p, re.I):
            out["order"] = p
        if re.search(r"viridae$", p, re.I):
            out["family"] = p
    last = parts[-1]
    if re.search(r"virus$", last, re.I) and not re.search(r"viridae$", last, re.I):
        out["genus"] = last
    return out


def load_meta() -> pd.DataFrame:
    def read_and_parse(df: pd.DataFrame):
        df.columns = df.columns.astype(str).str.replace("\ufeff", "", regex=False).str.strip()
        cols = {c.lower().strip(): c for c in df.columns}

        def pick(*names):
            for n in names:
                if n in cols:
                    return df[cols[n]]
            return None

        sid = pick("sample_id", "sampleid", "samplename", "run", "run_accession")
        if sid is None:
            return None, None
        out = pd.DataFrame({"sample_id": sid.astype(str).str.strip()})
        country = pick("country", "nation")
        if country is not None:
            out["Country"] = country.astype(str).str.strip()
        else:
            out["Country"] = "NA"
        study = pick("study", "bioproject", "project", "project_accession")
        if study is None:
            study = pick("project")
        if study is not None:
            out["Study"] = study.astype(str).str.strip()
        else:
            out["Study"] = "NA"
        out = out.drop_duplicates(subset=["sample_id"], keep="first")
        return out, None

    if META_XLSX.exists():
        try:
            dfx = pd.read_excel(META_XLSX, engine="openpyxl")
            out, _ = read_and_parse(dfx)
            if out is not None:
                return out
        except Exception:
            pass
    df = pd.read_csv(META_TSV, sep="\t")
    out, _ = read_and_parse(df)
    if out is None:
        raise ValueError("No meta information of ID（sample_id）；please check xlsx/tsv colname")
    return out


def continent_from_country(country: str) -> str:
    c = str(country).strip()
    asia = {"China", "Japan", "South Korea", "Thailand"}
    oceania = {"Australia", "New Zealand"}
    na = {"USA", "Canada"}
    sa = {"Brazil"}
    eu = {
        "UK", "France", "Germany", "Spain", "Denmark", "Norway", "Austria",
        "Ireland",
    }
    africa = {"Ghana", "Gabon"}
    if c in asia:
        return "Asia"
    if c in oceania:
        return "Oceania"
    if c in na:
        return "North America"
    if c in sa:
        return "South America"
    if c in africa:
        return "Africa"
    if c in eu:
        return "Europe"
    return "Other"


def read_matrix_header(path: Path) -> tuple[list[str], list[str]]:
    with path.open("r", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
    if not header or header[0] != "Contig":
        raise ValueError("First column should be Contig")
    sample_cols = header[1:]
    sample_ids = [c[: -len("_Count")] if c.endswith("_Count") else c for c in sample_cols]
    return header, sample_ids


def apply_matplotlib_pdf_editable_arial(base_size: float = 10.5) -> None:
    """
    Embed TrueType in PDF (fonttype 42); in Adobe this is usually editable text, not outlined Type 3.
    Prefer Arial; fall back to Helvetica / DejaVu if Arial is unavailable.
    """
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


def apply_nature_style() -> None:
    import matplotlib as mpl

    mpl.rcParams.update(
        {
            "figure.dpi": 120,
            "savefig.dpi": 300,
            "font.family": "sans-serif",
            "font.sans-serif": [
                "Arial",
                "Helvetica",
                "Arial Unicode MS",
                "DejaVu Sans",
                "Liberation Sans",
            ],
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "axes.titleweight": "normal",
            "axes.linewidth": 0.6,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.color": "#E5E5E5",
            "grid.linewidth": 0.4,
            "grid.linestyle": "-",
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def safe_unlink(path: Path) -> None:
    """Try to delete existing file to avoid overwrite failure when Excel has it open."""
    p = Path(path)
    for _ in range(8):
        try:
            if p.exists():
                p.unlink()
            return
        except PermissionError:
            time.sleep(0.35)


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    df.to_csv(tmp, sep="\t", index=False)
    safe_unlink(path)
    try:
        tmp.replace(path)
    except OSError:
        safe_unlink(path)
        tmp.replace(path)


def save_figure_dual(fig, stem: Path):
    """stem without suffix: write both PDF and PNG."""
    import matplotlib.pyplot as plt

    stem.parent.mkdir(parents=True, exist_ok=True)
    p_pdf = stem.with_suffix(".pdf")
    p_png = stem.with_suffix(".png")
    safe_unlink(p_pdf)
    safe_unlink(p_png)
    fig.savefig(p_pdf, bbox_inches="tight", format="pdf", facecolor="white", edgecolor="none")
    fig.savefig(p_png, bbox_inches="tight", format="png", facecolor="white", edgecolor="none", dpi=300)
    plt.close(fig)
    print("已写:", p_pdf)
    print("已写:", p_png)


def build_country_mean_long(
    accum_rank: dict[str, np.ndarray],
    meta_countries: pd.DataFrame,
    sample_ids: list[str],
    rank_key: str,
) -> pd.DataFrame:
    """accum_rank: taxon -> length n_samples vector of counts"""
    if not accum_rank:
        return pd.DataFrame(columns=["Country", rank_key, "mean_prop"])
    taxa = sorted(accum_rank.keys())
    mat = np.stack([accum_rank[t] for t in taxa], axis=0)
    col_sum = mat.sum(axis=0, keepdims=True)
    col_sum[col_sum == 0] = 1.0
    prop = mat / col_sum
    comp_samples = pd.DataFrame(prop.T, columns=taxa)
    comp_samples.insert(0, "sample_id", sample_ids)
    comp_long = comp_samples.merge(meta_countries, on="sample_id", how="left")
    countries = [c for c in comp_long["Country"].dropna().unique() if str(c) != "NA"]
    rows = []
    for cty in countries:
        sub = comp_long[comp_long["Country"] == cty]
        mean_prop = sub[taxa].mean(axis=0)
        rows.append(
            pd.DataFrame({"Country": cty, rank_key: taxa, "mean_prop": mean_prop.values})
        )
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def plot_stacked_taxa(
    comp_cty: pd.DataFrame,
    rank_key: str,
    countries_ord: list[str],
    top_n: int,
    title_suffix: str,
) -> tuple:
    import matplotlib.pyplot as plt

    col_name = rank_key
    top_ids = (
        comp_cty.groupby(col_name)["mean_prop"]
        .sum()
        .sort_values(ascending=False)
        .head(top_n)
        .index.tolist()
    )
    comp_plot = comp_cty[comp_cty[col_name].isin(top_ids)].copy()
    other = (
        comp_cty[~comp_cty[col_name].isin(top_ids)]
        .groupby("Country", as_index=False)["mean_prop"]
        .sum()
    )
    other[col_name] = "Other"
    comp_plot = pd.concat([comp_plot, other], ignore_index=True)

    pivot = comp_plot.pivot_table(
        index="Country", columns=col_name, values="mean_prop", fill_value=0.0
    )
    pivot = pivot.reindex([c for c in countries_ord if c in pivot.index])

    fig_w = min(12.0, max(7.0, 0.35 * len(pivot.index) + 4))
    fig_h = 4.8
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    bottom = np.zeros(len(pivot))
    n_col = len(pivot.columns)
    colors = [NATURE_PALETTE[i % len(NATURE_PALETTE)] for i in range(n_col)]
    for i, col in enumerate(pivot.columns):
        ax.bar(
            pivot.index,
            pivot[col].values,
            bottom=bottom,
            label=col,
            color=colors[i],
            width=0.72,
            edgecolor="white",
            linewidth=0.4,
        )
        bottom += pivot[col].values
    ax.set_ylabel("Mean relative abundance")
    en = RANK_LABELS_EN.get(rank_key, rank_key.title())
    ax.set_title(f"Viral {en} composition by country\n({title_suffix})")
    ax.set_xlabel("")
    ax.grid(axis="x", visible=False)
    ax.legend(
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        fontsize=6.5,
        ncol=1,
        borderaxespad=0,
        handlelength=1.2,
        handletextpad=0.5,
    )
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    fig.tight_layout()
    return fig, pivot


def main():
    import matplotlib.pyplot as plt

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    apply_nature_style()

    meta = load_meta()
    meta["Continent"] = meta["Country"].map(continent_from_country)

    header, sample_ids = read_matrix_header(MATRIX)
    n_samples = len(sample_ids)

    print("样本列数:", n_samples)
    print("加载病毒 taxonomy（多层级）…")
    tax_df = pd.read_csv(VIRUS_SUMMARY, sep="\t", usecols=["seq_name", "taxonomy"])
    votu_ranks: dict[str, dict[str, str]] = {}
    for r in tax_df.itertuples(index=False):
        votu_ranks[str(r.seq_name)] = parse_all_ranks(r.taxonomy)

    default_r = _default_ranks()
    accum: dict[str, dict[str, np.ndarray]] = {rk: {} for rk in RANK_KEYS}

    def add_to_rank(rk: str, label: str, row: np.ndarray):
        d = accum[rk]
        if label not in d:
            d[label] = np.zeros(n_samples, dtype=np.float64)
        d[label] += row

    observed = np.zeros(n_samples, dtype=np.int64)
    sum_n = np.zeros(n_samples, dtype=np.float64)
    sum_nlogn = np.zeros(n_samples, dtype=np.float64)

    print("流式扫描 Count 矩阵 …")
    chunk_i = 0
    for chunk in pd.read_csv(
        MATRIX,
        sep="\t",
        chunksize=CHUNKSIZE,
        dtype=np.int32,
        converters={"Contig": str},
        low_memory=False,
    ):
        contigs = chunk["Contig"].astype(str).to_numpy()
        data = chunk.drop(columns=["Contig"]).to_numpy(dtype=np.float64, copy=False)

        observed += (data > 0).sum(axis=0)
        sum_n += data.sum(axis=0)
        d = data.astype(np.float64, copy=False)
        with np.errstate(divide="ignore", invalid="ignore"):
            term = np.where(d > 0, d * np.log(d), 0.0)
        sum_nlogn += term.sum(axis=0)

        for r, row in enumerate(data):
            ranks = votu_ranks.get(contigs[r], default_r)
            for rk in RANK_KEYS:
                add_to_rank(rk, ranks[rk], row)

        chunk_i += 1
        if chunk_i % 100 == 0:
            print(f"  … {chunk_i * CHUNKSIZE} 行")

    shannon = np.zeros(n_samples, dtype=np.float64)
    mask = sum_n > 0
    shannon[mask] = np.log(sum_n[mask]) - sum_nlogn[mask] / sum_n[mask]

    alpha = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "observed_vOTUs": observed.astype(int),
            "shannon_natural": shannon,
            "viral_read_counts_sum": sum_n,
        }
    )
    alpha = alpha.merge(meta, on="sample_id", how="left")
    write_tsv(alpha, OUT_DIR / "biogeo_alpha_per_sample.tsv")
    print("已写:", OUT_DIR / "biogeo_alpha_per_sample.tsv")

    def grp_summary(df, key):
        return (
            df.groupby(key, dropna=False)
            .agg(
                n_samples=("sample_id", "count"),
                observed_median=("observed_vOTUs", "median"),
                observed_mean=("observed_vOTUs", "mean"),
                shannon_median=("shannon_natural", "median"),
                shannon_mean=("shannon_natural", "mean"),
                viral_counts_sum_median=("viral_read_counts_sum", "median"),
            )
            .reset_index()
        )

    for key, name in [
        ("Country", "biogeo_alpha_summary_by_country.tsv"),
        ("Continent", "biogeo_alpha_summary_by_continent.tsv"),
        ("Study", "biogeo_alpha_summary_by_study.tsv"),
    ]:
        sub = alpha[alpha[key].notna() & (alpha[key].astype(str) != "NA")]
        if len(sub):
            write_tsv(grp_summary(sub, key), OUT_DIR / name)
            print("已写:", OUT_DIR / name)

    cov = (
        alpha.groupby("Country", dropna=False)
        .agg(n_samples=("sample_id", "count"))
        .reset_index()
        .sort_values("n_samples", ascending=False)
    )
    write_tsv(cov, OUT_DIR / "biogeo_sampling_intensity_by_country.tsv")
    countries_ord = cov["Country"].tolist()

    meta_sub = meta[["sample_id", "Country"]].copy()

    # Country-mean composition tables and stacked charts per taxonomic level
    for rk in RANK_KEYS:
        long_df = build_country_mean_long(accum[rk], meta_sub, sample_ids, rk)
        tsv_name = f"biogeo_composition_{rk}_country_mean.tsv"
        write_tsv(long_df, OUT_DIR / tsv_name)
        print("已写:", OUT_DIR / tsv_name)

        if long_df.empty:
            continue
        top_n = TOP_N_BY_RANK.get(rk, 12)
        fig, _ = plot_stacked_taxa(
            long_df,
            rk,
            countries_ord,
            top_n,
            "mean of per-sample relative abundances; rare taxa as Other",
        )
        stem = OUT_DIR / f"biogeo_fig_{rk}_stacked_by_country"
        save_figure_dual(fig, stem)

    # —— Sampling intensity (styled) ——
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    order = cov["Country"].tolist()
    counts = cov.set_index("Country")["n_samples"].reindex(order).values
    x = np.arange(len(order))
    ax.bar(x, counts, color="#4E79A7", width=0.72, edgecolor="white", linewidth=0.4)
    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=45, ha="right")
    ax.set_ylabel("Number of samples")
    ax.set_title("Sampling intensity by country")
    ax.grid(axis="x", visible=False)
    fig.tight_layout()
    save_figure_dual(fig, OUT_DIR / "biogeo_fig_sampling_by_country")

    # —— Alpha diversity boxplots (styled) ——
    plot_alpha = alpha[
        alpha["Country"].notna() & (alpha["Country"].astype(str) != "NA")
    ].copy()
    plot_alpha["Country"] = pd.Categorical(
        plot_alpha["Country"], categories=countries_ord, ordered=True
    )

    cty_colors = country_color_map(countries_ord)

    fig, axes = plt.subplots(1, 2, figsize=(7.8, 3.8))
    for ax, col, ylab in zip(
        axes,
        ["observed_vOTUs", "shannon_natural"],
        ["Observed vOTUs (richness)", "Shannon diversity (ln)"],
    ):
        labels: list[str] = []
        data_by: list[np.ndarray] = []
        for c in countries_ord:
            d = plot_alpha.loc[plot_alpha["Country"] == c, col].values
            if len(d) > 0:
                labels.append(str(c))
                data_by.append(d)
        bp = ax.boxplot(
            data_by,
            labels=labels,
            patch_artist=True,
            widths=0.65,
            medianprops={"color": "#333333", "linewidth": 1.0},
            boxprops={"linewidth": 0.6},
            whiskerprops={"linewidth": 0.6},
            capprops={"linewidth": 0.6},
            flierprops={"marker": "o", "markersize": 2, "alpha": 0.35},
        )
        for i, cty in enumerate(labels):
            col_hex = cty_colors.get(cty, "#4E79A7")
            fill = lighten_hex_color(col_hex)
            bp["boxes"][i].set_facecolor(fill)
            bp["boxes"][i].set_edgecolor(col_hex)
            bp["medians"][i].set_color("#333333")
            for j in (2 * i, 2 * i + 1):
                bp["whiskers"][j].set_color(col_hex)
                bp["caps"][j].set_color(col_hex)
        ax.set_ylabel(ylab)
        ax.set_xlabel("")
        ax.grid(axis="x", visible=False)
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    fig.tight_layout()
    save_figure_dual(fig, OUT_DIR / "biogeo_fig_alpha_by_country")

    # Legacy phylum filename is included in the loop; loop already writes biogeo_fig_phylum_stacked_by_country
    print("Finished")


if __name__ == "__main__":
    main()
