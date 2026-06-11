# GSVC paper figure scripts

Standalone plotting scripts for the Global Swine Viral Catalogue manuscript. Each file is named `Fig{panel}_{short_description}.{R|py}` and writes PDFs (and often PNGs) under **`output/`**.

## Setup

Set the project root (folder containing `FINAL_all_projects_Count_matrix.tsv`):

```bash
# Linux / macOS
export PHAGE_PROJECT_ROOT="/path/to/Phage"

# Windows PowerShell
$env:PHAGE_PROJECT_ROOT = "D:\F\微生物Meta\global\Figure\Phage"
```

Run from anywhere:

```bash
Rscript paper_figure_scripts/Fig5g_lab_prev_ratio_boxplot.R
python paper_figure_scripts/Fig5a_prevalence_hist_ecdf.py
```

## Regenerating copies from upstream sources

Most scripts are adapted copies of the analysis repository:

```bash
python paper_figure_scripts/generate_paper_scripts.py
```

See `GENERATION_MANIFEST.tsv` for source → destination mapping.

## Figure index (main text)

| Panel | Script | Primary output | Notes |
|-------|--------|----------------|-------|
| **Fig. 1a** | `Fig1a_country_sample_map.R` | world map PDF | Needs `China_sheng.geojson` in working dir or edit path in script |
| **Fig. 1c** | `Fig1c_vOTU_length_ridge.R` | `output/Fig1c_vOTU_length_distribution.pdf` | Uses `high_medium_quality_length.txt` |
| **Fig. 1d** | `Fig1d_novelty_ICTV_IMGVR.R` | ICTV/IMG/VR bars | Embedded summary data |
| **Fig. 1e** | `Fig1e_database_comparison.R` | database comparison | |
| **Fig. 1f** | `Fig1f_read_capture.R` | read-capture rates | |
| **Fig. 1g** | `Fig1g_novel_vOTU_bars.R` | novel vOTU counts | |
| **Fig. 2a** | `Fig2a_taxonomy_composition_GSVC.R` | taxonomy panels a–b | Args: `--tsv`, `--out_dir` |
| **Fig. 2b** | `Fig2b_multi_family_circular_phylo.R` | circular multi-family tree | Requires vContact3 inputs |
| **Fig. 3a** | `Fig3a_alpha_by_country.py` | `03_Biogeographical/biogeo_fig_alpha_by_country.*` | Full biogeo descriptive module |
| **Fig. 3b** | `Fig3b_rarefaction_by_country.py` | `biogeo_fig_rarefaction_by_country.*` | |
| **Fig. 3c** | `Fig3c_pcoa_bray_country_centroids.py` | PCoA centroids | Large PERMANOVA script |
| **Fig. 3d/e** | `Fig3d_distance_decay_within_country.py` | decay / hexbin Mantel | |
| **Fig. 3f** | `Fig3f_neutral_model_sloan.py` | Sloan-like neutral model | Pass `--outdir` |
| **Fig. 4a–b** | `Fig4a_host_phylum_genus_bars.R` | host taxonomy top25 | |
| **Fig. 4c–d** | `Fig4c_host_range_breadth_lines.R` | host range / infection load | Same source file |
| **Fig. 4e** | `Fig4e_abundance_vOTU_correlation.R` | abundance vs vOTU richness | |
| **Fig. 4g** | `Fig4g_kill_the_winner_stratified.R` | stratified phage tiers | |
| **Fig. 4h** | `Fig4h_age_binscatter_species.R` | age-stratified binscatter | |
| **Fig. 5a–b** | `Fig5a_prevalence_hist_ecdf.py` | prevalence hist + ECDF | |
| **Fig. 5c** | `Fig5c_core_tpm_by_gut_stacked.R` | core TPM by gut | `--plot-only` optional |
| **Fig. 5d** | `Fig5d_ge05_host_genus_by_gut.R` | Lacto/Prevotella-like by gut | Run `analyze_ge05_votu_host_genus_tpm_by_gut.py` first if TSV missing |
| **Fig. 5e** | `Fig5e_PCoA_PAM_k4.R` or `Fig5e_PCoA_PAM_k4_plot_only.R` | `PCoA_PAM_k4.pdf` | Full vs copy existing |
| **Fig. 5f** | `Fig5f_enterotype_driver_heatmap.R` | driver heatmap | Needs PAM assignments |
| **Fig. 5g** | `Fig5g_lab_prev_ratio_boxplot.R` | Lab/Prev log2 ratio | |
| **Fig. 5h** | `Fig5h_enterotype_age_pies.R` | age pies by enterotype | |
| **Fig. 6a** | `Fig6a_eggnog_top_terms.R` | eggNOG top terms | See Functional outputs |
| **Fig. 6b** | `Fig6b_ARG_burden_by_country.R` | ARG burden by country | |
| **Fig. 6c** | `Fig6c_AMG_host_metabolism_heatmap.R` | AMG × host heatmap | |
| **Fig. 6d–e** | `Fig6d_jaccard_null_model.R`, `Fig6e_KO_directionality_forest.py` | KO coupling | Run host_phage_ko pipeline first |
| **Fig. 7a** | `Fig7a_treatment_schematic.py` | cohort summary | |
| **Fig. 7b** | `Fig7b_phage_alpha_lineplot_7arm.R` | alpha diversity lines | PRJNA1010706 |
| **Fig. 7c** | `Fig7c_ARG_phage_fraction_facets.R` | ARG phage fraction | |
| **Fig. 7d** | `Fig7d_lsalivarius_prophage_scatter.py` | L. salivarius scatter | |
| **Fig. 7e** | `Fig7e_provirus_linear_map.R` | provirus linear map | |

## Supplementary figures (selection)

| Panel | Script |
|-------|--------|
| Fig. S2 | `FigS2_vOTU_length_study_compare.R` |
| Fig. S3 | `FigS3_cross_study_clustering.R` |
| Fig. S5 | `FigS5_lifestyle_pie.R` |
| Fig. S6–S8 | `FigS6_KO_single_vs_broad_host.R`, `FigS7_CAZy_single_vs_broad.R`, `FigS8_integrase_prevalence.R` |
| Fig. S9 | `FigS9_lifestyle_abundance_correlation.R` |
| Fig. S10 | `FigS10_core_taxonomy_pies.py` |
| Fig. S11 | `FigS11_enterotype_alpha_boxplot.R` |
| Fig. S12–S15 | `FigS12_AMG_metabolism_bars.R` … `FigS15_AMG_distribution.R` |
| Fig. S17–S22 | `FigS17_shannon_by_arm_facets.R` … `FigS22_lsalivarius_evidence_chain.py` |

## Enterotype reproduction (Fig. 5e–h)

Publication parameters used in the manuscript:

```bash
cd /path/to/Phage
Rscript paper_figure_scripts/Fig5e_PCoA_PAM_k4.R \
  --feces-only --min-depth 10000 --prev 0.30 --top-taxa 5000 \
  --k-max 4 --main-grp NONE --out-dir Enterotype_vOTU_k4_2
```

Then:

```bash
Rscript paper_figure_scripts/Fig5f_enterotype_driver_heatmap.R --assign meta_qc/Enterotype_vOTU_k4_2/PAM_k4_samples.tsv ...
Rscript paper_figure_scripts/Fig5g_lab_prev_ratio_boxplot.R --assign meta_qc/PAM_k4_samples.tsv ...
```

## Not yet included (add manually if needed)

- **Fig. 1a / 1b**: metadata overview panels (custom ggplot from `meta_augmented_combined.tsv`)
- **Fig. S1**: CheckV QC summary (derive from `high_medium_quality_checkv_results.tsv`)
- **Fig. S4**: geNomad length × taxonomy (`Fig_geNomad_class_and_family_stacked_by_length` pipeline)
- **Fig. 4f**: host range × lifestyle bar (extract section from `host_range_and_infection_stats.R`)

## Dependencies

- **R**: ggplot2, data.table, vegan, cluster, philentropy, ggpubr, ggridges, patchwork, circlize (figure-dependent)
- **Python**: pandas, numpy, scipy, matplotlib, statsmodels

Versions are listed in the manuscript Key Resources Table.
