# GSVC_github publication script
#!/usr/bin/env python3
"""Copy and adapt repository figure scripts into paper_figure_scripts/Fig*.R|py."""
from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "_common"))
from encoding_utils import read_text_auto, write_text_utf8  # noqa: E402

PHAGE = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent

R_HEADER = r'''# =============================================================================
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

'''

PY_HEADER = '''# =============================================================================
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

'''

# (dest_name, relative_source_from_PHAGE, panel_note)
FIGURES = [
    # Figure 1 (finalFigures/Figure1)
    ("Fig1a_country_sample_map.R", "finalFigures/Figure1/country_sample_counts_by_continent_positions.R", "Fig. 1a world map"),
    ("Fig1c_vOTU_length_ridge.R", "finalFigures/Figure1/Fig 1c.R", "Fig. 1c length ridgeline"),
    ("Fig1d_novelty_ICTV_IMGVR.R", "finalFigures/Figure1/Fig 1d.R", "Fig. 1d novelty bars"),
    ("Fig1e_database_comparison.R", "finalFigures/Figure1/Fig 1e.R", "Fig. 1e database comparison"),
    ("Fig1f_read_capture.R", "finalFigures/Figure1/Fig 1f.R", "Fig. 1f read capture"),
    ("Fig1g_novel_vOTU_bars.R", "finalFigures/Figure1/Fig 1g.R", "Fig. 1g novel vOTUs"),
    ("FigS2_vOTU_length_study_compare.R", "finalFigures/Figure1/Fig S2/Distribution of vOTU Lengths(study对比）/优化.R", "Fig. S2"),
    ("FigS3_cross_study_clustering.R", "finalFigures/Figure1/Fig S3.R", "Fig. S3"),
    # Figure 2
    ("Fig2a_taxonomy_composition_GSVC.R", "Functional/scripts/plot_fig3_taxonomy_composition_R.R", "Fig. 2a"),
    # Figure 3
    ("Fig3a_alpha_by_country.py", "biogeography_descriptive_analysis.py", "Fig. 3a (runs full biogeo module)"),
    ("Fig3b_rarefaction_by_country.py", "03_Biogeographical/biogeography_rarefaction_by_country.py", "Fig. 3b"),
    ("Fig3c_pcoa_bray_country_centroids.py", "04_CLR_PERMANOVA/compositional_continent_analysis.py", "Fig. 3c"),
    ("Fig3d_distance_decay_within_country.py", "04_CLR_PERMANOVA/distance_decay_spatial_subsets.py", "Fig. 3d/e"),
    ("Fig3f_neutral_model_sloan.py", "make_paper_figs_part2.py", "Fig. 3f"),
    # Figure 4
    ("Fig4a_host_phylum_genus_bars.R", "host_taxonomy_vOTU_top20_plot.R", "Fig. 4a–b"),
    ("Fig4c_host_range_breadth_lines.R", "host_range_and_infection_stats.R", "Fig. 4c–d"),
    ("Fig4e_abundance_vOTU_correlation.R", "abundance_lifestyle_vOTU_correlation.R", "Fig. 4e"),
    ("Fig4f_host_range_lifestyle_bar.R", "host_range_and_infection_stats.R", "Fig. 4f (same source; see README)"),
    ("Fig4g_kill_the_winner_stratified.R", "abundance_Virulent_vOTU_correlation.R", "Fig. 4g"),
    ("Fig4h_age_binscatter_species.R", "binscatter_overlay_stage_two_species_feces.R", "Fig. 4h"),
    # Figure 5
    ("Fig5a_prevalence_hist_ecdf.py", "06_core_virome_figs/core_prevalence_taxonomy_host_pies.py", "Fig. 5a–b"),
    ("Fig5c_core_tpm_by_gut_stacked.R", "meta_qc/plot_core_votu_tpm_fraction_by_gut.R", "Fig. 5c"),
    ("Fig5d_ge05_host_genus_by_gut.R", "meta_qc/plot_ge05_votu_host_genus_by_gut_stacked.R", "Fig. 5d"),
    ("Fig5e_PCoA_PAM_k4.R", "meta_qc/run_enterotype_votu.R", "Fig. 5e (full enterotype pipeline)"),
    ("Fig5f_enterotype_driver_heatmap.R", "meta_qc/plot_enterotype_marker_heatmap.R", "Fig. 5f"),
    ("Fig5g_lab_prev_ratio_boxplot.R", "meta_qc/plot_enterotype_lab_prev_ratio_boxplot.R", "Fig. 5g"),
    ("Fig5h_enterotype_age_pies.R", "meta_qc/plot_enterotype_cluster_age_pie.R", "Fig. 5h"),
    ("FigS10_core_taxonomy_pies.py", "06_core_virome_figs/core_prevalence_taxonomy_host_pies.py", "Fig. S10"),
    ("FigS11_enterotype_alpha_boxplot.R", "meta_qc/plot_enterotype_alpha_diversity.R", "Fig. S11"),
    # Figure 6
    ("Fig6a_eggnog_top_terms.R", "Functional/plot_vibrant_amg_metabolism_distribution.R", "Fig. 6a (see README)"),
    ("Fig6b_ARG_burden_by_country.R", "Functional/plot_ARG_country_strict_vs_nonstrict.R", "Fig. 6b"),
    ("Fig6c_AMG_host_metabolism_heatmap.R", "Functional/plot_vibrant_metabolism_host_heatmap_FDR.R", "Fig. 6c"),
    ("Fig6d_jaccard_null_model.R", "Functional/host_phage_ko_coupling_figures.py", "Fig. 6d"),
    ("Fig6e_KO_directionality_forest.py", "Functional/host_phage_ko_coupling_figures.py", "Fig. 6e"),
    ("FigS12_AMG_metabolism_bars.R", "Functional/plot_vibrant_amg_metabolism_distribution.R", "Fig. S12"),
    ("FigS13_ARG_nonstrict_country.R", "Functional/plot_ARG_country_strict_vs_nonstrict.R", "Fig. S13"),
    ("FigS14_VF_box_by_country.R", "Functional/vf_phage_host_plots.R", "Fig. S14"),
    ("FigS15_AMG_distribution.R", "Functional/plot_vibrant_amg_metabolism_distribution.R", "Fig. S15"),
    ("FigS5_lifestyle_pie.R", "lifestyle_virulent_temperate_proportion.R", "Fig. S5"),
    ("FigS6_KO_single_vs_broad_host.R", "meta_qc/PRJNA1010706_treatment_analysis/plot_KO_enrichment_single_vs_broad_bar.R", "Fig. S6"),
    ("FigS7_CAZy_single_vs_broad.R", "meta_qc/PRJNA1010706_treatment_analysis/plot_host_range_CAZy_integrase_bar.R", "Fig. S7"),
    ("FigS8_integrase_prevalence.R", "meta_qc/PRJNA1010706_treatment_analysis/plot_host_range_CAZy_integrase_bar.R", "Fig. S8"),
    ("FigS9_lifestyle_abundance_correlation.R", "abundance_lifestyle_vOTU_correlation.R", "Fig. S9"),
    # Figure 7
    ("Fig2b_multi_family_circular_phylo.R", "Functional/scripts/plot_multi_family_circular_phylo.R", "Fig. 2b"),
    ("Fig7a_treatment_schematic.py", "meta_qc/PRJNA1010706_meta_merged/plot_meta_summary.py", "Fig. 7a cohort schematic"),
    ("Fig7b_phage_alpha_lineplot_7arm.R", "meta_qc/PRJNA1010706_treatment_analysis/plot_lineplot_7arm_4age_quick.R", "Fig. 7b"),
    ("Fig7c_ARG_phage_fraction_facets.R", "meta_qc/PRJNA1010706_treatment_analysis/plot_arg_frac_boxline_facets.R", "Fig. 7c"),
    ("Fig7d_lsalivarius_prophage_scatter.py", "meta_qc/PRJNA1010706_treatment_analysis/26_plot_lsalivarius_scatter_pub.py", "Fig. 7d"),
    ("Fig7e_provirus_linear_map.R", "meta_qc/PRJNA1010706_treatment_analysis/plot_provirus_gggenes_SRR15732359.R", "Fig. 7e"),
    ("FigS17_shannon_by_arm_facets.R", "meta_qc/PRJNA1010706_treatment_analysis/plot_shannon_boxline_facets.R", "Fig. S17"),
    ("FigS18_pcoa_treatment_age.R", "meta_qc/PRJNA1010706_treatment_analysis/plot_pcoa_treatment_age_quick.R", "Fig. S18"),
    ("FigS19_ARG_frac_vs_G6.R", "meta_qc/PRJNA1010706_treatment_analysis/plot_arg_frac_vs_G6_by_age_wilcox.R", "Fig. S19"),
    ("FigS20_KO_heatmap_treatment.R", "meta_qc/PRJNA1010706_treatment_analysis/plot_functional_kegg_cazy_heatmap.R", "Fig. S20"),
    ("FigS21_CAZy_heatmap_treatment.R", "meta_qc/PRJNA1010706_treatment_analysis/plot_functional_kegg_cazy_heatmap.R", "Fig. S21"),
    ("FigS22_lsalivarius_evidence_chain.py", "meta_qc/PRJNA1010706_treatment_analysis/25_lsalivarius_prophage_evidence_no_mediation.py", "Fig. S22"),
]

REPLACEMENTS = [
    (r'[Ff]:[/\\]微生物Meta[/\\]global[/\\]Figure[/\\]Phage', 'PHAGE_ROOT'),
    (r'[Ff]:[/\\]微生物Meta[/\\]global[/\\]Figure[/\\]Figure1', 'file.path(PHAGE_ROOT, "finalFigures", "Figure1")'),
    (r'[Dd]:[/\\]F[/\\]微生物Meta[/\\]global[/\\]Figure[/\\]Phage', 'PHAGE_ROOT'),
    (r'setwd\s*\(\s*"[^"]+"\s*\)', '# setwd: using PHAGE_ROOT above'),
    (r'if \(dir\.exists\("D:\\\\F\\\\微生物Meta\\\\global\\\\Figure\\\\Phage"\)\) setwd\("D:\\\\F\\\\微生物Meta\\\\global\\\\Figure\\\\Phage"\)',
     '# setwd: using PHAGE_ROOT above'),
]


def adapt_body(text: str, is_r: bool) -> str:
    for pat, rep in REPLACEMENTS:
        if is_r and rep == "PHAGE_ROOT":
            text = re.sub(pat, '"', text)  # skip broken replace
            continue
        text = re.sub(pat, rep, text)
    # Fix R PHAGE_ROOT string paths
    text = text.replace('"PHAGE_ROOT"', 'PHAGE_ROOT')
    text = re.sub(
        r'F:[/\\]微生物Meta[/\\]global[/\\]Figure[/\\]Phage',
        'PHAGE_ROOT',
        text,
    )
    text = re.sub(
        r'D:[/\\]F[/\\]微生物Meta[/\\]global[/\\]Figure[/\\]Phage',
        'PHAGE_ROOT',
        text,
    )
    text = re.sub(
        r'F:[/\\]微生物Meta[/\\]global[/\\]Figure[/\\]Figure1',
        'file.path(PHAGE_ROOT, "finalFigures", "Figure1")',
        text,
    )
    return text


def main() -> None:
    manifest = []
    for dest, src_rel, note in FIGURES:
        src = PHAGE / src_rel
        dest_path = OUT / dest
        if not src.exists():
            manifest.append(f"MISSING\t{dest}\t{src_rel}\t{note}")
            continue
        body = read_text_auto(src)
        is_r = dest.endswith(".R")
        body = adapt_body(body, is_r)
        header = R_HEADER if is_r else PY_HEADER
        write_text_utf8(dest_path, header + "\n" + body)
        manifest.append(f"OK\t{dest}\t{src_rel}\t{note}")

    (OUT / "GENERATION_MANIFEST.tsv").write_text(
        "status\tdest\tsource\tnote\n" + "\n".join(manifest) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {sum(1 for m in manifest if m.startswith('OK'))} scripts; "
          f"{sum(1 for m in manifest if m.startswith('MISSING'))} missing.")


if __name__ == "__main__":
    main()
