# Script index (main figures)

| Figure | Script directory | Key scripts |
|--------|------------------|-------------|
| Fig. 1 | `Fig01_catalogue/` | `Fig1a`–`Fig1g`, `FigS2`, `FigS3` |
| Fig. 2 | `Fig02_taxonomy/` | `Fig2a`, `Fig2b` |
| Fig. 3 | `Fig03_biogeography/` | `Fig3a`–`Fig3f` |
| Fig. 4 | `Fig04_host_phage/` | `Fig4a`–`Fig4h` |
| Fig. 5 | `Fig05_core_enterotype/` | `Fig5a`–`Fig5h`, `FigS10`, `FigS11` |
| Fig. 6 | `Fig06_function/` | `Fig6a`–`Fig6e`, `FigS12`–`FigS15` |
| Fig. 7 | `Fig07_intervention/` | `Fig7a`–`Fig7e`, `FigS17`–`FigS22` |

## Quick start

```bash
export GSVC_ROOT=/path/to/GSVC_github
cd $GSVC_ROOT

# Example: Fig. 4a host phylum bars
Rscript Scripts/Fig04_host_phage/Fig4a_host_phylum_genus_bars.R

# Example: Fig. 5a prevalence (requires Count matrix from Zenodo 10.5281/zenodo.20579313)
python Scripts/Fig05_core_enterotype/Fig5a_prevalence_hist_ecdf.py

# Example: Fig. 7b (redraw from cohort tables)
cd meta_qc/PRJNA1010706_treatment_analysis
Rscript ../../Scripts/Fig07_intervention/Fig7b_phage_alpha_lineplot_7arm.R .
```

Outputs go to `Scripts/output/` or cohort `figures/` subfolders.

See `README.md` in this folder for the full panel-level table.
