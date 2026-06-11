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
# vOTU enterotype analysis: pvca_Y_count.tsv + meta_augmented_combined.tsv
#
# Core: relative abundance → zero-safe JSD distance → PAM + hclust (default)
# Optional: DMM (DirichletMultinomial; may be slow on high-dimensional sparse data)
#
# Example run (Phage root):
#   Rscript meta_qc/run_enterotype_votu.R
#   Rscript meta_qc/run_enterotype_votu.R --feces-only --main-grp Age_group --min-depth 3000 --prev 0.10 --top-taxa 3000
#   Rscript meta_qc/run_enterotype_votu.R --max-samples 1500   # trial: subsample for distance
#   Rscript meta_qc/run_enterotype_votu.R --run-dmm            # optional: run DMM
#   Rscript meta_qc/run_enterotype_votu.R --no-marginal        # skip marginal density
#   Rscript meta_qc/run_enterotype_votu.R --biplot-n 10        # vOTU loading arrows if needed
#   Rscript meta_qc/run_enterotype_votu.R --marginal-size 12 # flatter marginals (ggMarginal size)
#
# Output dir: meta_qc/Enterotype_vOTU/ (--out-dir to change)
# After main run, meta_qc/plot_followup_enterotype_run.R with same count/prev/depth
# can plot alpha diversity + marker heatmap (followup_figures/ under out-dir).
#
# -----------------------------------------------------------------------------
# Parameters vs cluster separation (no need to grid-search full data blindly)
#
# 1) Bottleneck is sample count n: JSD is O(n^2). Tune with --max-samples 2000–4000 first
#    (--feces-only) for quick PCoA + pam_k_silhouette.tsv; then full run without max-samples.
#
# 2) --prev (min prevalence): higher → fewer vOTUs, distance on high-prevalence core, less noise,
#    but low-prevalence drivers may drop. 0.50 is aggressive; try 0.05–0.20.
#
# 3) --top-taxa: after prev, keep top N by mean rel. abundance; too large → sparse noise, slow;
#    too small → lose signal. Often tune with prev (e.g. 5000–15000 when prev is low).
#
# 4) --min-depth: higher → fewer, deeper-sequenced samples; cleaner silhouette, selection bias.
#
# 5) --k-max / --pam-k: silhouette picks k. Use --pam-k 3 to fix k=3,
#    skip searching 2–10. Use --k-max 4 for quick screening.
#
# 6) --main-grp: if column missing, no shape; color by enterotype only.
#
# Batch tuning: enterotype_param_screen.example.sh (subsample + grid to out-dirs).
# -----------------------------------------------------------------------------
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

parse_args <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  getv <- function(flag, default = NA_character_) {
    i <- match(flag, a)
    if (is.na(i) || i >= length(a)) return(default)
    a[[i + 1]]
  }
  has <- function(flag) flag %in% a
  list(
    out_dir = getv("--out-dir", "Enterotype_vOTU"),
    main_grp = getv("--main-grp", "Age_group"),
    sample_id = getv("--sample-id", "sample_id"),
    count_path = getv("--count-path", ""),
    meta_path = getv("--meta-path", ""),
    gut_col = getv("--gut-col", "Gut_harmonized"),
    gut_keep = getv("--gut-keep", "feces"),
    feces_only = has("--feces-only"),
    min_depth = as.integer(getv("--min-depth", "3000")),
    prev = as.numeric(getv("--prev", "0.05")),
    top_taxa = as.integer(getv("--top-taxa", "3000")),
    max_samples = {
      s <- getv("--max-samples", "")
      if (!nzchar(s)) NA_integer_ else as.integer(s)
    },
    seed = as.integer(getv("--seed", "123")),
    k_min = as.integer(getv("--k-min", "2")),
    k_max = as.integer(getv("--k-max", "10")),
    pam_k = {
      s <- getv("--pam-k", "")
      if (!nzchar(s)) NA_integer_ else as.integer(s)
    },
    run_dmm = has("--run-dmm"),
    eps = as.numeric(getv("--eps", "1e-6")),
    marginal = !has("--no-marginal"),
    biplot_n = as.integer(getv("--biplot-n", "0")),
    # ggMarginal: size = main vs marginal ratio; larger → shorter marginals (CRAN default 5)
    marginal_size = as.integer(getv("--marginal-size", "10"))
  )
}

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
root <- if (length(script_path)) dirname(dirname(normalizePath(script_path))) else getwd()
if (!file.exists(file.path(root, "pvca_Y_count.tsv"))) root <- getwd()

pa <- parse_args()

count_path <- if (nzchar(pa$count_path)) pa$count_path else file.path(root, "pvca_Y_count.tsv")
meta_path  <- if (nzchar(pa$meta_path)) pa$meta_path else file.path(root, "meta_qc", "meta_augmented_combined.tsv")
if (!file.exists(count_path)) stop("No: ", count_path)
if (!file.exists(meta_path)) stop("No: ", meta_path)

out_dir <- file.path(root, "meta_qc", pa$out_dir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

meta <- fread(meta_path, sep = "\t", encoding = "UTF-8", na.strings = c("", "NA", "na", "NaN"))
meta[, (pa$sample_id) := as.character(get(pa$sample_id))]

# Feces-only samples (recommended for enterotype)
if (isTRUE(pa$feces_only)) {
  g <- tolower(trimws(as.character(meta[[pa$gut_col]])))
  keep <- g == tolower(pa$gut_keep)
  keep[is.na(keep)] <- FALSE
  meta <- meta[keep]
}

cnt <- fread(count_path, sep = "\t", encoding = "UTF-8", showProgress = TRUE)
cnt[, sample_id := as.character(sample_id)]
if (anyDuplicated(cnt$sample_id)) {
  cnt <- unique(cnt, by = "sample_id")
}

common <- intersect(meta[[pa$sample_id]], cnt$sample_id)
meta <- meta[get(pa$sample_id) %in% common]
cnt  <- cnt[sample_id %in% common]

set.seed(pa$seed)
if (is.finite(pa$max_samples) && !is.na(pa$max_samples) && pa$max_samples > 0L && nrow(meta) > pa$max_samples) {
  keep_ids <- sample(meta[[pa$sample_id]], pa$max_samples)
  meta <- meta[get(pa$sample_id) %in% keep_ids]
  cnt  <- cnt[sample_id %in% keep_ids]
}

# Align row order
setkey(cnt, sample_id)
cnt <- cnt[J(meta[[pa$sample_id]])]
stopifnot(all(cnt$sample_id == meta[[pa$sample_id]]))

feat_cols <- setdiff(names(cnt), "sample_id")
mat <- as.matrix(cnt[, ..feat_cols])
storage.mode(mat) <- "double"
rownames(mat) <- cnt$sample_id
colnames(mat) <- feat_cols

# Sample depth filtering
depth <- rowSums(mat)
keep_s <- depth >= pa$min_depth
mat <- mat[keep_s, , drop = FALSE]
meta <- meta[match(rownames(mat), meta[[pa$sample_id]])]

# Feature prevalence filtering
prev <- colSums(mat > 0) / nrow(mat)
keep_f <- prev >= pa$prev
mat <- mat[, keep_f, drop = FALSE]

# Take top taxa (by mean relative abundance)
rel0 <- mat / rowSums(mat)
mean_rel <- colMeans(rel0)
if (is.finite(pa$top_taxa) && !is.na(pa$top_taxa) && pa$top_taxa > 0L && ncol(mat) > pa$top_taxa) {
  ord <- order(mean_rel, decreasing = TRUE)
  take <- ord[seq_len(pa$top_taxa)]
  mat <- mat[, take, drop = FALSE]
  rel0 <- rel0[, take, drop = FALSE]
  mean_rel <- colMeans(rel0)
}


# JSD: add eps pseudocount and normalize to avoid Inf/NaN from zeros
eps <- pa$eps
rel <- rel0 + eps
rel <- rel / rowSums(rel)

dist_mat <- as.matrix(philentropy::distance(rel, method = "jensen-shannon", unit = "log2"))
sid <- rownames(rel)
rownames(dist_mat) <- sid
colnames(dist_mat) <- sid
dist_jsd <- as.dist(dist_mat)
saveRDS(dist_jsd, file.path(out_dir, "dist_jsd.rds"))

# PCoA
pcoa <- cmdscale(dist_jsd, k = 3, eig = TRUE)
pcoa_df <- as.data.frame(pcoa$points)
names(pcoa_df) <- c("PCo1", "PCo2", "PCo3")
pcoa_df$sample_id <- rownames(pcoa_df)
var_exp <- round(pcoa$eig / sum(pcoa$eig) * 100, 1)

meta_df <- as.data.frame(meta)
meta_df$sample_id <- meta_df[[pa$sample_id]]
pcoa_df <- merge(pcoa_df, meta_df, by = "sample_id", all.x = TRUE, sort = FALSE)

# PAM: choose best k (silhouette)
k_min <- pa$k_min
k_max <- pa$k_max
ks <- k_min:k_max

pam_eval <- data.frame(k = ks, silhouette = NA_real_)
if (is.finite(pa$pam_k) && !is.na(pa$pam_k) && pa$pam_k >= 2L) {
  best_k <- pa$pam_k
} else {
  for (k in ks) {
    pfit <- cluster::pam(dist_jsd, k = k, diss = TRUE)
    sil <- cluster::silhouette(pfit$clustering, dist_jsd)
    pam_eval$silhouette[pam_eval$k == k] <- mean(sil[, "sil_width"])
  }
  best_k <- pam_eval$k[which.max(pam_eval$silhouette)]
}
fwrite(as.data.table(pam_eval), file.path(out_dir, "pam_k_silhouette.tsv"), sep = "\t")

pam_fit <- cluster::pam(dist_jsd, k = best_k, diss = TRUE)
pam_assign <- data.frame(sample_id = names(pam_fit$clustering), cluster = as.integer(pam_fit$clustering))
fwrite(as.data.table(pam_assign), file.path(out_dir, sprintf("PAM_k%d_samples.tsv", best_k)), sep = "\t")

# hclust (k from PAM best k)
hc <- hclust(dist_jsd, method = "ward.D2")
hc_assign <- data.frame(sample_id = hc$labels, cluster = as.integer(cutree(hc, k = best_k)))
fwrite(as.data.table(hc_assign), file.path(out_dir, sprintf("HCLUST_k%d_samples.tsv", best_k)), sep = "\t")
saveRDS(hc, file.path(out_dir, "hclust_wardD2.rds"))

# Driver vOTUs (PAM): top 10 mean rel. abundance per cluster
assign2 <- pam_assign[match(rownames(rel0), pam_assign$sample_id), , drop = FALSE]
if (anyNA(assign2$sample_id)) {
  missing_n <- sum(is.na(assign2$sample_id))
}
keep_idx <- !is.na(assign2$sample_id)
rel0_k <- rel0[keep_idx, , drop = FALSE]
assign2 <- assign2[keep_idx, , drop = FALSE]
drivers <- rbindlist(lapply(sort(unique(assign2$cluster)), function(k) {
  m <- colMeans(rel0_k[assign2$cluster == k, , drop = FALSE])
  o <- order(m, decreasing = TRUE)
  topn <- seq_len(min(10L, length(o)))
  data.table(cluster = k, taxon = names(m)[o][topn], mean_rel = m[o][topn])
}))
fwrite(drivers, file.path(out_dir, sprintf("PAM_k%d_driver_votu_top10.tsv", best_k)), sep = "\t")

# PCoA plot (color by cluster; shape by main-grp; optional marginal density + vOTU arrows)
suppressPackageStartupMessages({
  library(ggplot2)
})
plot_df <- merge(pcoa_df, pam_assign, by = "sample_id", all.x = TRUE, sort = FALSE)
plot_df <- plot_df[!is.na(plot_df$cluster), , drop = FALSE]

main_grp <- pa$main_grp
use_shape <- main_grp %in% names(plot_df)
if (isTRUE(use_shape)) {
  plot_df[[main_grp]] <- as.factor(plot_df[[main_grp]])
  plot_df <- plot_df[!is.na(plot_df[[main_grp]]), , drop = FALSE]
} else {
  message("meta no include --main-grp column ", main_grp)
}

# —— Biplot: Spearman of vOTU rel. abundance vs PCo1/2; top N arrows by |vector| (exploratory) ——
arrows_df <- NULL
bn <- pa$biplot_n
if (is.finite(bn) && !is.na(bn) && bn > 0L && nrow(plot_df) >= 5L) {
  rel_plot <- rel0[match(plot_df$sample_id, rownames(rel0)), , drop = FALSE]
  if (nrow(rel_plot) == nrow(plot_df) && ncol(rel_plot) > 0L) {
    c1 <- suppressWarnings(apply(rel_plot, 2L, function(v) {
      cor(v, plot_df$PCo1, method = "spearman", use = "complete.obs")
    }))
    c2 <- suppressWarnings(apply(rel_plot, 2L, function(v) {
      cor(v, plot_df$PCo2, method = "spearman", use = "complete.obs")
    }))
    c1[is.na(c1)] <- 0
    c2[is.na(c2)] <- 0
    magn <- sqrt(c1^2 + c2^2)
    take <- order(magn, decreasing = TRUE)[seq_len(min(bn, length(magn)))]
    take <- take[magn[take] > 0]
    if (length(take) > 0L) {
      span <- max(diff(range(plot_df$PCo1)), diff(range(plot_df$PCo2)), na.rm = TRUE)
      sc <- 0.42 * span / max(magn[take], 1e-9)
      arrows_df <- data.frame(
        taxon = names(magn)[take],
        cor_PCo1 = c1[take],
        cor_PCo2 = c2[take],
        magn = magn[take],
        x1 = c1[take] * sc,
        y1 = c2[take] * sc,
        stringsAsFactors = FALSE
      )
      fwrite(as.data.table(arrows_df), file.path(out_dir, sprintf("PCoA_PAM_k%d_biplot_arrows.tsv", best_k)), sep = "\t")
    }
  }
}

# Enterotype colors: matplotlib tab20 (match host genus stacked bars)
cl_ids <- sort(unique(as.integer(plot_df$cluster)))
TAB20_COL <- c(
  "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a",
  "#d62728", "#ff9896", "#9467bd", "#c5b0d5", "#8c564b", "#c49c94",
  "#e377c2", "#f7b6d2", "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d",
  "#17becf", "#9edae5"
)
pal_ent <- TAB20_COL[((seq_along(cl_ids) - 1L) %% length(TAB20_COL)) + 1L]
names(pal_ent) <- as.character(cl_ids)

if (isTRUE(use_shape)) {
  p <- ggplot(plot_df, aes(PCo1, PCo2, color = factor(cluster), shape = .data[[main_grp]]))
} else {
  p <- ggplot(plot_df, aes(PCo1, PCo2, color = factor(cluster)))
}

# Layer order: ellipse, arrows, points (no vOTU text labels)
p <- p +
  stat_ellipse(
    data = plot_df,
    aes(PCo1, PCo2, group = factor(cluster)),
    inherit.aes = FALSE,
    linetype = 2,
    linewidth = 0.4
  )

if (!is.null(arrows_df) && nrow(arrows_df) > 0L) {
  p <- p +
    geom_segment(
      data = arrows_df,
      aes(x = 0, y = 0, xend = x1, yend = y1),
      inherit.aes = FALSE,
      linewidth = 0.35,
      color = "grey35",
      arrow = grid::arrow(length = grid::unit(0.12, "cm"), type = "closed")
    )
}

if (isTRUE(use_shape)) {
  p <- p +
    geom_point(size = 2.35, alpha = 0.68) +
    guides(shape = guide_legend(override.aes = list(alpha = 1)))
} else {
  p <- p + geom_point(size = 2.35, alpha = 0.68, stroke = 0)
}

p <- p +
  scale_color_manual(
    values = pal_ent,
    breaks = as.character(cl_ids),
    name = "Enterotype",
    drop = FALSE
  ) +
  coord_fixed() +
  theme_classic(base_size = 12) +
  theme(
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10),
    legend.key.height = unit(0.6, "lines"),
    axis.title = element_text(size = 11),
    axis.text  = element_text(size = 10),
    plot.title = element_text(size = 12, hjust = 0.5),
    plot.margin = margin(5.5, 12, 5.5, 5.5)
  ) +
  labs(
    color = "Enterotype",
    x = paste0("PCo1 (", var_exp[1], "%)"),
    y = paste0("PCo2 (", var_exp[2], "%)"),
    title = sprintf("vOTU enterotype (PAM, k=%d)", best_k)
  )

if (isTRUE(use_shape)) {
  p <- p + labs(shape = main_grp)
}

# Marginal density (color by enterotype); requires ggExtra
w_pdf <- 7.8
h_pdf <- 6.0
ms <- pa$marginal_size
if (is.na(ms) || ms < 2L) ms <- 5L

if (isTRUE(pa$marginal) && requireNamespace("ggExtra", quietly = TRUE)) {
  pm <- ggExtra::ggMarginal(
    p,
    margins = "both",
    type = "density",
    size = ms,
    groupColour = TRUE,
    groupFill = TRUE
  )
  ggplot2::ggsave(file.path(out_dir, sprintf("PCoA_PAM_k%d.pdf", best_k)), pm, width = w_pdf, height = h_pdf, useDingbats = FALSE)
  ggplot2::ggsave(file.path(out_dir, sprintf("PCoA_PAM_k%d.png", best_k)), pm, width = w_pdf, height = h_pdf, dpi = 300)
} else {
  if (isTRUE(pa$marginal)) {
  }
  ggsave(file.path(out_dir, sprintf("PCoA_PAM_k%d.pdf", best_k)), p, width = 7.2, height = 5.6, useDingbats = FALSE)
  ggsave(file.path(out_dir, sprintf("PCoA_PAM_k%d.png", best_k)), p, width = 7.2, height = 5.6, dpi = 300)
}

# Save version without marginals (layout / legacy)
ggsave(file.path(out_dir, sprintf("PCoA_PAM_k%d_no_marginal.pdf", best_k)), p, width = 7.2, height = 5.6, useDingbats = FALSE)
ggsave(file.path(out_dir, sprintf("PCoA_PAM_k%d_no_marginal.png", best_k)), p, width = 7.2, height = 5.6, dpi = 300)

  # DirichletMultinomial typically uses sample×feature count matrix
  count_mat <- mat
  k_max <- pa$k_max
  lplc <- rep(NA_real_, k_max)
  fit_list <- vector("list", k_max)
  for (k in seq_len(k_max)) {
    fit_list[[k]] <- DirichletMultinomial::dmn(count = count_mat, k = k, verbose = FALSE)
    lplc[[k]] <- DirichletMultinomial::laplace(fit_list[[k]])
  }
  dmm_eval <- data.frame(k = seq_len(k_max), laplace = lplc)
  fwrite(as.data.table(dmm_eval), file.path(out_dir, "dmm_k_laplace.tsv"), sep = "\t")
  best_k_dmm <- which.min(lplc)
  z <- DirichletMultinomial::mixture(fit_list[[best_k_dmm]])
  dmm_assign <- data.frame(sample_id = rownames(z), cluster = apply(z, 1, which.max))
  fwrite(as.data.table(dmm_assign), file.path(out_dir, sprintf("DMM_k%d_samples.tsv", best_k_dmm)), sep = "\t")
}

