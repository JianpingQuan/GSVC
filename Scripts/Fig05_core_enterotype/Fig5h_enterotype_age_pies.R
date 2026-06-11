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
# Age distribution within enterotype — pie charts
#
# Two age sources (--age-bins):
#   group  use existing meta category (default Age_group), see --age-col
#   days   use meta Age_days_parsed, fixed day bins:
#          Lactation 0–21d | Nursery 22–63d | Growing 64–119d |
#          Early finishing 120–180d | Late finishing >180d; invalid/missing → unlabeled
#
# Scope:
#   Default: one --cluster (no sector labels; legend optional)
#   --all-clusters: 2×2 panel (four ETs), shared legend, no pie labels;
#                   still write _ETk_counts.tsv and _all_clusters_counts.tsv
#
# Usage (all four ETs + five age bins + combined):
#   Rscript meta_qc/plot_enterotype_cluster_age_pie.R \
#     --assign meta_qc/PAM_k4_samples.tsv \
#     --all-clusters \
#     --age-bins days \
#     --meta meta_qc/meta_augmented_combined.tsv \
#     --out-prefix meta_qc/enterotype_k4_downloaded/enterotype_age_days_pies
#
#   Outputs: {out_prefix}_four_panel.pdf/.png, _ET*_counts.tsv, _all_clusters_counts.tsv
#
# Requires: data.table, ggplot2
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

DAY_LEVELS <- c(
  "Lactation (0–21d)",
  "Nursery (22–63d)",
  "Growing (64–119d)",
  "Early finishing (120–180d)",
  "Late finishing (>180d)",
  "未标注 / 缺失"
)

# Match DAY_LEVELS order (shared legend colors)
DAY_COLORS <- c(
  "#66c2a5",
  "#fc8d62",
  "#8da0cb",
  "#e78ac3",
  "#a6d854",
  "#999999"
)

bin_age_days <- function(x) {
  d <- suppressWarnings(as.numeric(x))
  lab <- rep(NA_character_, length(d))
  bad <- !is.finite(d) | d < 0
  lab[bad] <- DAY_LEVELS[[6L]]
  ok <- !bad
  lab[ok & d <= 21] <- DAY_LEVELS[[1L]]
  lab[ok & d >= 22 & d <= 63] <- DAY_LEVELS[[2L]]
  lab[ok & d >= 64 & d <= 119] <- DAY_LEVELS[[3L]]
  lab[ok & d >= 120 & d <= 180] <- DAY_LEVELS[[4L]]
  lab[ok & d > 180] <- DAY_LEVELS[[5L]]
  still <- ok & is.na(lab)
  if (any(still)) lab[still] <- DAY_LEVELS[[6L]]
  factor(lab, levels = DAY_LEVELS)
}

parse_args <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  getv <- function(flag, default = NA_character_) {
    i <- match(flag, a)
    if (is.na(i) || i >= length(a)) return(default)
    a[[i + 1]]
  }
  list(
    assign_path = getv("--assign", ""),
    meta_path = getv("--meta", ""),
    cluster = as.integer(getv("--cluster", "3")),
    age_col = getv("--age-col", "Age_group"),
    age_bins = tolower(trimws(getv("--age-bins", "group"))),
    out_prefix = getv("--out-prefix", ""),
    all_clusters = "--all-clusters" %in% a
  )
}

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
root <- if (length(script_path)) dirname(dirname(normalizePath(script_path))) else getwd()
if (!file.exists(file.path(root, "pvca_Y_count.tsv"))) root <- getwd()

pa <- parse_args()


assign_path <- if (file.exists(pa$assign_path)) normalizePath(pa$assign_path) else file.path(root, pa$assign_path)

meta_path <- if (nzchar(pa$meta_path)) {
  if (file.exists(pa$meta_path)) normalizePath(pa$meta_path) else file.path(root, pa$meta_path)
} else {
  file.path(root, "meta_qc", "meta_augmented_combined.tsv")
}


out_prefix <- if (nzchar(pa$out_prefix)) {
  normalizePath(pa$out_prefix, winslash = "/", mustWork = FALSE)
} else {
  file.path(dirname(assign_path), "enterotype_age_pie")
}
dir.create(dirname(out_prefix), showWarnings = FALSE, recursive = TRUE)

ac <- pa$age_col

asg <- fread(assign_path, sep = "\t", encoding = "UTF-8")

asg[, cluster := as.integer(cluster)]
asg[, sample_id := trimws(as.character(sample_id))]


meta <- fread(meta_path, sep = "\t", encoding = "UTF-8")
meta[, sample_id := trimws(as.character(sample_id))]

if (pa$age_bins == "days") {
  if (!"Age_days_parsed" %in% names(meta)) {
    stop()
  }
  meta_sub <- unique(meta[, .(sample_id, age_days = Age_days_parsed)], by = "sample_id")
  merged_all <- merge(asg, meta_sub, by = "sample_id", all.x = TRUE)
  merged_all[, age_lab := bin_age_days(age_days)]
  age_title_bit <- "Age（d）"
  legend_title <- "Age"
} else {
  if (!ac %in% names(meta)) {
    stop()
  }
  meta_sub <- unique(meta[, .(sample_id, age_val = get(ac))], by = "sample_id")
  merged_all <- merge(asg, meta_sub, by = "sample_id", all.x = TRUE)
  merged_all[, age_lab := {
    v <- trimws(as.character(age_val))
    v[!nzchar(v) | v %in% c("NA", "N/A", "na", "n/a")] <- "No labels"
    factor(v)
  }]
  age_title_bit <- paste0("Classfy：", ac)
  legend_title <- ac
}

clusters <- if (isTRUE(pa$all_clusters)) {
  sort(unique(asg$cluster))
} else {
  pa$cluster
}

assign_bn <- basename(assign_path)
meta_bn <- basename(meta_path)
long_rows <- list()
panel_rows <- list()

pal_base <- c(
  "#8dd3c7", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69",
  "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5"
)

# Single-ET pie (no sector labels); day mode uses shared colors + legend
plot_one_pie <- function(tab, k, n_merged, age_title_bit, assign_bn, meta_bn, show_legend) {
  tab <- copy(tab)
  tab[, age_lab := as.character(age_lab)]
  tab[, pct := 100 * N / sum(N)]
  if (pa$age_bins == "days") {
    tab[, age_lab := factor(age_lab, levels = DAY_LEVELS)]
    tab <- tab[order(age_lab)]
    fill_map <- setNames(DAY_COLORS, DAY_LEVELS)
    fill_vals <- fill_map[as.character(tab$age_lab)]
  } else {
    tab <- tab[order(-N)]
    ulev <- unique(tab$age_lab)
    fill_vals <- pal_base[seq_along(ulev)]
    names(fill_vals) <- as.character(ulev)
  }

  ggplot(tab, aes(x = "", y = N, fill = age_lab)) +
    geom_col(width = 1, color = "white", linewidth = 0.35) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = fill_vals, drop = FALSE, name = legend_title) +
    labs(
      title = paste0("ET", k, "：", age_title_bit),
      subtitle = paste0("n = ", n_merged, " ；", assign_bn),
      x = NULL,
      y = NULL,
      caption = meta_bn
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 8.5, colour = "grey30"),
      plot.caption = element_text(hjust = 0.5, size = 7.5, colour = "grey40"),
      legend.position = if (isTRUE(show_legend)) "bottom" else "none",
      legend.title = element_text(face = "bold", size = 9),
      plot.margin = margin(10, 10, 10, 10)
    )
}

for (k in clusters) {
  merged <- merged_all[cluster == k]
  if (!nrow(merged)) {
    message()
    next
  }
  tab <- merged[, .N, by = age_lab]
  stem <- if (isTRUE(pa$all_clusters)) paste0(out_prefix, "_ET", k) else out_prefix
  fwrite(tab[, .(age_lab, N, pct = 100 * N / sum(N))], paste0(stem, "_counts.tsv"), sep = "\t")
  long_rows[[length(long_rows) + 1L]] <- tab[, .(cluster = k, age_lab, N)]

  if (!isTRUE(pa$all_clusters)) {
    p <- plot_one_pie(tab, k, nrow(merged), age_title_bit, assign_bn, meta_bn, show_legend = TRUE)
    ggsave(paste0(stem, ".pdf"), p, width = 7.8, height = 7.4, limitsize = FALSE)
    ggsave(paste0(stem, ".png"), p, width = 7.8, height = 7.4, dpi = 300, limitsize = FALSE)
  } else {
    tr <- copy(tab)
    if (pa$age_bins == "days") {
      tr[, age_lab := factor(age_lab, levels = DAY_LEVELS)]
    }
    tr[, et := factor(paste0("ET", k), levels = paste0("ET", clusters))]
    panel_rows[[length(panel_rows) + 1L]] <- tr[, .(et, age_lab, N)]
  }
}

if (length(long_rows)) {
  allc <- rbindlist(long_rows)
  allc[, pct := 100 * N / sum(N), by = cluster]
  fwrite(allc, paste0(out_prefix, "_all_clusters_counts.tsv"), sep = "\t")
}

if (isTRUE(pa$all_clusters) && length(panel_rows)) {
  comb <- rbindlist(panel_rows)
  et_levels <- paste0("ET", clusters)
  comb[, et := factor(et, levels = et_levels)]
  if (pa$age_bins == "days") {
    comb[, age_lab := factor(age_lab, levels = DAY_LEVELS)]
    fill_scale <- scale_fill_manual(
      values = setNames(DAY_COLORS, DAY_LEVELS),
      breaks = DAY_LEVELS,
      drop = FALSE,
      name = legend_title
    )
  } else {
    u <- sort(unique(as.character(comb$age_lab)))
    comb[, age_lab := factor(age_lab, levels = u)]
    fill_scale <- scale_fill_manual(
      values = setNames(colorRampPalette(pal_base)(length(u)), u),
      limits = u,
      drop = FALSE,
      name = legend_title
    )
  }

  p_panel <- ggplot(comb, aes(x = "", y = N, fill = age_lab)) +
    geom_col(width = 1, color = "white", linewidth = 0.3) +
    coord_polar(theta = "y") +
    facet_wrap(~et, nrow = 2L, ncol = 2L, scales = "free_y") +
    fill_scale +
    labs(
      subtitle = paste0("assign: ", assign_bn, ""),
      caption = paste0("meta: ", meta_bn),
      x = NULL,
      y = NULL
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "grey32"),
      plot.caption = element_text(hjust = 0.5, size = 8, colour = "grey40"),
      strip.text = element_text(face = "bold", size = 11),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 8.5),
      legend.key.height = unit(0.55, "cm"),
      plot.margin = margin(10, 12, 10, 10)
    )

  ns <- merged_all[, .(n = .N), keyby = cluster][, paste(paste0("ET", cluster, "=", n), collapse = "  ")]
  p_panel <- p_panel + labs(subtitle = paste0(assign_bn, "  |  ", ns))

  ggsave(
    paste0(out_prefix, "_four_panel.pdf"),
    p_panel,
    width = 10.2,
    height = 8.6,
    limitsize = FALSE
  )
  ggsave(
    paste0(out_prefix, "_four_panel.png"),
    p_panel,
    width = 10.2,
    height = 8.6,
    dpi = 300,
    limitsize = FALSE
  )
}
