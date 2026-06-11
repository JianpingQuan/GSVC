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


#!/usr/bin/env Rscript
# Quick redraw: phage_votu_alpha_diversity_lineplot_7arm_4age (from saved tables).
suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
analysis_dir <- if (length(args) >= 1L) {
  normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

source(file.path(analysis_dir, "R", "science_aaas_figures.R"))
theme_set(theme_science_aaas(11))

res_dir <- file.path(analysis_dir, "results", "phage_votu_diversity")
fig_dir <- file.path(analysis_dir, "figures", "phage_votu_diversity")
age_levels <- c(19L, 23L, 34L, 48L)
arm_levels_7 <- paper_arm_levels()
arm_cols_7 <- paper_arm_palette()

summ <- read_tsv(file.path(res_dir, "alpha_diversity_mean_se_by_arm_age.tsv"), show_col_types = FALSE) |>
  mutate(
    paper_arm = factor(paper_arm, levels = arm_levels_7),
    age_d = factor(as.integer(as.character(age_d)), levels = age_levels),
    metric = factor(metric, levels = c("Shannon", "Observed richness"))
  )

kw <- read_tsv(file.path(res_dir, "kruskal_wallis_alpha_7arm_by_age.tsv"), show_col_types = FALSE) |>
  mutate(
    age_d = factor(as.integer(age_d), levels = age_levels),
    metric = factor(metric, levels = c("Shannon", "Observed richness"))
  )

p_to_star <- function(p) {
  if (length(p) != 1L || !is.finite(p)) {
    return(NA_character_)
  }
  if (p < 0.001) {
    return("***")
  }
  if (p < 0.01) {
    return("**")
  }
  if (p < 0.05) {
    return("*")
  }
  NA_character_
}

ann_y <- summ |>
  group_by(age_d, metric) |>
  summarise(y_top = max(mean + se, na.rm = TRUE), .groups = "drop")

ann_stars <- kw |>
  mutate(stars = vapply(kw_p, p_to_star, character(1))) |>
  filter(!is.na(stars)) |>
  left_join(ann_y, by = c("age_d", "metric")) |>
  mutate(y = y_top * 1.1)

p_line <- ggplot(summ, aes(x = age_d, y = mean, colour = paper_arm, group = paper_arm)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.14, linewidth = 0.45) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2.2) +
  geom_text(
    data = ann_stars,
    aes(x = age_d, y = y, label = stars),
    inherit.aes = FALSE,
    size = 4.6,
    fontface = "bold",
    colour = "#1A1A1A",
    vjust = 0
  ) +
  facet_wrap(~metric, scales = "free_y", ncol = 1L) +
  scale_colour_manual(values = arm_cols_7, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.22))) +
  scale_x_discrete(labels = function(x) paste0(x, " d")) +
  labs(
    x = "Sampling age (d)",
    y = "Alpha diversity",
    title = "Phage vOTU alpha diversity over time (seven treatment arms)",
    subtitle = "GG + G1\u2013G6; mean \u00b1 SE; stars = Kruskal\u2013Wallis across arms at each age",
    caption = "Stars: * p<0.05, ** p<0.01, *** p<0.001 (seven-group Kruskal\u2013Wallis, nominal)"
  ) +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

ggsave(
  file.path(fig_dir, "phage_votu_alpha_diversity_lineplot_7arm_4age.pdf"),
  p_line,
  width = 4.8,
  height = 8.4,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(fig_dir, "phage_votu_alpha_diversity_lineplot_7arm_4age.png"),
  p_line,
  width = 4.8,
  height = 8.4,
  dpi = 320,
  bg = "white"
)

message("Wrote lineplot: ", fig_dir)
