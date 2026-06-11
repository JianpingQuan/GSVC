# GSVC_github publication script
# =============================================================================
# Fig. 5e — PCoA of fecal vOTU enterotypes (PAM k = 4), plot-only
# Requires prior enterotype run OR bundled coordinates in meta_qc/.
# Full pipeline: paper_figure_scripts/Fig5e_PCoA_PAM_k4.R
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

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# Search enterotype output directories (newest first)
cand_dirs <- c(
  file.path(PHAGE_ROOT, "meta_qc", "Enterotype_vOTU_k4_2"),
  file.path(PHAGE_ROOT, "meta_qc", "Enterotype_vOTU"),
  file.path(PHAGE_ROOT, "meta_qc")
)
k <- 4L
pcoa_tsv <- NA_character_
pam_tsv <- file.path(PHAGE_ROOT, "meta_qc", sprintf("PAM_k%d_samples.tsv", k))
for (d in cand_dirs) {
  f <- file.path(d, sprintf("PCoA_PAM_k%d_biplot_arrows.tsv", k))
  if (file.exists(file.path(d, sprintf("PCoA_PAM_k%d.pdf", k)))) {
    pcoa_pdf <- file.path(d, sprintf("PCoA_PAM_k%d.pdf", k))
    file.copy(pcoa_pdf, file.path(OUT_DIR, "Fig5e_PCoA_PAM_k4.pdf"), overwrite = TRUE)
    message("Copied existing PCoA: ", pcoa_pdf)
    quit(save = "no", status = 0)
  }
}

# Re-run enterotype with publication parameters (slow; omit for quick plot if PDF exists above)
message("No precomputed PCoA PDF found. Run:\n",
        "  Rscript paper_figure_scripts/Fig5e_PCoA_PAM_k4.R \\\n",
        "    --feces-only --min-depth 10000 --prev 0.30 --top-taxa 5000 --k-max 4 \\\n",
        "    --out-dir Enterotype_vOTU_k4_2")
stop("Fig5e plot-only: missing PCoA_PAM_k4.pdf in meta_qc enterotype folders.")
