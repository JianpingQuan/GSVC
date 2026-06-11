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


# =============================================================================
# Redraw Fig_core_votu_tpm_fraction_by_gut_stacked (R / ggplot2)
# Same logic as analyze_core_votu_tpm_by_gut.py; legend top horizontal, larger font.
#
# Usage (project root with FINAL_all_projects_TPM_matrix.tsv):
#   Rscript meta_qc/plot_core_votu_tpm_fraction_by_gut.R
#   Rscript meta_qc/plot_core_votu_tpm_fraction_by_gut.R --plot-only   # plot from existing core_votu_tpm_fraction_by_gut.tsv
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
root <- if (nzchar(script_path)) {
  dirname(dirname(normalizePath(script_path)))
} else {
  getwd()
}
if (!file.exists(file.path(root, "FINAL_all_projects_TPM_matrix.tsv"))) {
  root <- getwd()
}

pa <- commandArgs(trailingOnly = TRUE)
plot_only <- "--plot-only" %in% pa

OUT_DIR <- file.path(root, "meta_qc")
TPM_MATRIX <- file.path(root, "FINAL_all_projects_TPM_matrix.tsv")
METRICS_PATH <- file.path(root, "core_virome_vOTU_metrics.csv")
META_PATH <- file.path(root, "meta_qc", "meta_augmented_combined.tsv")
TSV_WIDE <- file.path(OUT_DIR, "core_votu_tpm_fraction_by_gut.tsv")

GUT_ORDER <- c("Jejunum", "Ileum", "cecum", "Colon", "feces")
PREV_STRICT <- 0.9
PREV_LIKE_LO <- 0.5
CHUNK_ROWS <- 800L

gut_xlabel <- function(g) {
  s <- trimws(as.character(g))
  out <- s
  len <- nchar(s)
  sub1 <- len == 1L & nzchar(s)
  sub2 <- len > 1L
  out[sub1] <- toupper(s[sub1])
  ii <- which(sub2)
  if (length(ii)) {
    out[ii] <- paste0(toupper(substr(s[ii], 1L, 1L)), substr(s[ii], 2L, nchar(s[ii])))
  }
  out
}

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

if (isTRUE(plot_only)) {
  if (!file.exists(TSV_WIDE)) stop("找不到 ", TSV_WIDE, "；请先跑 Python 分析或去掉 --plot-only")
  wide <- as.data.frame(fread(TSV_WIDE, sep = "\t"))
  ord <- match(GUT_ORDER, wide$Gut_harmonized)
  wide <- wide[ord[!is.na(ord)], , drop = FALSE]
} else {
  if (!file.exists(TPM_MATRIX)) stop("找不到: ", TPM_MATRIX)
  if (!file.exists(METRICS_PATH)) stop("找不到: ", METRICS_PATH)
  if (!file.exists(META_PATH)) stop("找不到: ", META_PATH)

  meta <- fread(META_PATH, sep = "\t", encoding = "UTF-8", na.strings = c("", "NA"))
  meta <- meta[!is.na(sample_id) & !is.na(Gut_harmonized)]
  meta[, sample_id := trimws(as.character(sample_id))]
  meta[, Gut_harmonized := trimws(as.character(Gut_harmonized))]
  meta <- unique(meta, by = "sample_id")
  meta <- meta[Gut_harmonized != "" & Gut_harmonized != "Cecal mucosal" & Gut_harmonized != "Stomach"]
  id_to_gut <- setNames(meta$Gut_harmonized, meta$sample_id)

  hdr <- fread(TPM_MATRIX, sep = "\t", nrows = 0L)
  raw_cols <- names(hdr)
  first_col <- raw_cols[[1L]]
  stems <- vapply(raw_cols[-1L], function(cn) {
    cn <- as.character(cn)
    if (endsWith(cn, "_TPM")) substr(cn, 1L, nchar(cn) - 4L) else cn
  }, character(1L))
  keep <- stems[stems %in% names(id_to_gut)]
  if (length(keep) < 1L) stop("TPM 列与 meta sample_id 无交集")
  tpm_col_names <- paste0(keep, "_TPM")
  gut_per_col <- unname(id_to_gut[keep])
  n_s <- length(keep)

  prev_dt <- fread(METRICS_PATH, sep = ",", select = c("vOTU", "prevalence"))
  prev_dt[, `:=`(vOTU = trimws(as.character(vOTU)), prevalence = as.numeric(prevalence))]
  prev_dt <- unique(prev_dt, by = "vOTU")

  strict_sum <- numeric(n_s)
  like_sum <- numeric(n_s)
  other_sum <- numeric(n_s)

  col_names <- c(first_col, tpm_col_names)
  offset <- 0L
  repeat {
    dt <- fread(
      TPM_MATRIX,
      sep = "\t",
      skip = 1L + offset,
      nrows = CHUNK_ROWS,
      header = FALSE,
      col.names = col_names,
      showProgress = TRUE
    )
    if (nrow(dt) == 0L) break
    Contig <- trimws(as.character(dt[[1L]]))
    mat <- as.matrix(dt[, -1L, with = FALSE])
    storage.mode(mat) <- "double"
    pv <- prev_dt$prevalence[match(Contig, prev_dt$vOTU)]
    pv[is.na(pv)] <- 0
    ge_s <- pv >= PREV_STRICT
    ge_l <- pv >= PREV_LIKE_LO & pv < PREV_STRICT
    ge_o <- pv < PREV_LIKE_LO
    strict_sum <- strict_sum + colSums(sweep(mat, 1L, ge_s, `*`))
    like_sum <- like_sum + colSums(sweep(mat, 1L, ge_l, `*`))
    other_sum <- other_sum + colSums(sweep(mat, 1L, ge_o, `*`))
    offset <- offset + nrow(dt)
  }

  total <- strict_sum + like_sum + other_sum
  pct_s <- ifelse(total > 0, strict_sum / total * 100, NA_real_)
  pct_l <- ifelse(total > 0, like_sum / total * 100, NA_real_)
  pct_o <- ifelse(total > 0, other_sum / total * 100, NA_real_)

  rows <- list()
  for (g in GUT_ORDER) {
    ix <- which(gut_per_col == g)
    if (length(ix) == 0L) next
    rows[[length(rows) + 1L]] <- data.frame(
      Gut_harmonized = g,
      mean_pct_TPM_core_strict = mean(pct_s[ix], na.rm = TRUE),
      mean_pct_TPM_core_like = mean(pct_l[ix], na.rm = TRUE),
      mean_pct_TPM_non_core = mean(pct_o[ix], na.rm = TRUE),
      n_samples = length(ix),
      stringsAsFactors = FALSE
    )
  }
  wide <- do.call(rbind, rows)

  fwrite(as.data.table(wide), TSV_WIDE, sep = "\t")
  long_dt <- rbindlist(lapply(seq_len(nrow(wide)), function(i) {
    r <- wide[i, ]
    data.table(
      Gut_harmonized = r$Gut_harmonized,
      category = c(
        "Core-strict (prev \u2265 0.9)",
        "Core-like (0.5 \u2264 prev < 0.9)",
        "Non-core (prev < 0.5)"
      ),
      mean_pct_TPM = c(r$mean_pct_TPM_core_strict, r$mean_pct_TPM_core_like, r$mean_pct_TPM_non_core),
      n_samples = r$n_samples
    )
  }))
  fwrite(long_dt, file.path(OUT_DIR, "core_votu_tpm_fraction_by_gut_long.tsv"), sep = "\t")
  fwrite(wide[, c("Gut_harmonized", "n_samples")], file.path(OUT_DIR, "core_votu_tpm_fraction_by_gut_sample_counts.tsv"), sep = "\t")
}

# ---- ggplot: stacked bars, legend top horizontal, large font ----
gut_levels <- gut_xlabel(wide$Gut_harmonized)
plot_long <- rbind(
  data.frame(
    Gut = factor(gut_levels, levels = gut_levels),
    Category = "Core-strict (prev \u2265 0.9)",
    pct = wide$mean_pct_TPM_core_strict,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Gut = factor(gut_levels, levels = gut_levels),
    Category = "Core-like (0.5 \u2264 prev < 0.9)",
    pct = wide$mean_pct_TPM_core_like,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Gut = factor(gut_levels, levels = gut_levels),
    Category = "Non-core (prev < 0.5)",
    pct = wide$mean_pct_TPM_non_core,
    stringsAsFactors = FALSE
  )
)
plot_long$Category <- factor(
  plot_long$Category,
  levels = c("Core-strict (prev \u2265 0.9)", "Core-like (0.5 \u2264 prev < 0.9)", "Non-core (prev < 0.5)")
)

fills <- c(
  "Core-strict (prev \u2265 0.9)" = "#2E6F95",
  "Core-like (0.5 \u2264 prev < 0.9)" = "#8FBC8F",
  "Non-core (prev < 0.5)" = "#D4D4D4"
)

p <- ggplot(plot_long, aes(x = Gut, y = pct, fill = Category)) +
  geom_col(width = 0.62, color = "white", linewidth = 0.35) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_fill_manual(values = fills, name = NULL) +
  labs(x = NULL, y = "Mean fraction of vOTU TPM (%)") +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.box = "horizontal",
    plot.title = element_blank(),
    axis.text = element_text(size = 15, colour = "black"),
    axis.title.y = element_text(size = 17, margin = margin(r = 8)),
    legend.text = element_text(size = 15),
    legend.margin = margin(b = 4),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(12, 14, 10, 10)
  ) +
  guides(fill = guide_legend(nrow = 1L, byrow = TRUE))

ggsave(file.path(OUT_DIR, "Fig_core_votu_tpm_fraction_by_gut_stacked.pdf"), p, width = 9.5, height = 6.2, useDingbats = FALSE)
ggsave(file.path(OUT_DIR, "Fig_core_votu_tpm_fraction_by_gut_stacked.png"), p, width = 9.5, height = 6.2, dpi = 300)

message("Finished: ", normalizePath(file.path(OUT_DIR, "Fig_core_votu_tpm_fraction_by_gut_stacked.pdf"), winslash = "/", mustWork = FALSE))
