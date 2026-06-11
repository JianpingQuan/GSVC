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
# =============================================================================
# ARG phage fraction: faceted boxplots by paper_arm; colour = cohort_Age.
# Friedman (4 ages, complete pig blocks) from pipeline TSV.
#
# Usage:
#   Rscript plot_arg_frac_boxline_facets.R
#   Rscript plot_arg_frac_boxline_facets.R "D:/path/to/PRJNA1010706_treatment_analysis"
#
# Requires: ggplot2, readr, dplyr
# =============================================================================

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

arm_levels <- c("GG", "G1", "G2", "G3", "G4", "G5", "G6")
age_levels <- c(19L, 23L, 34L, 48L)

source(file.path(analysis_dir, "R", "science_aaas_figures.R"), local = TRUE)

res_dir <- file.path(analysis_dir, "results", "arg_phage_sarg_fraction")
fig_dir <- file.path(analysis_dir, "figures", "arg_phage_sarg_fraction")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

samp_path <- file.path(res_dir, "sample_arg_phage_fraction.tsv")
fr_path <- file.path(res_dir, "friedman_time_within_arm.tsv")
if (!file.exists(samp_path)) {
  stop("Missing ", samp_path, " — run 16_arg_phage_fraction_treatment.py first.")
}

df <- read_tsv(samp_path, show_col_types = FALSE) |>
  filter(is.finite(arg_frac)) |>
  mutate(
    cohort_Age = as.integer(round(as.numeric(cohort_Age))),
    Age_f = factor(cohort_Age, levels = age_levels),
    paper_arm = factor(paper_arm, levels = arm_levels)
  ) |>
  filter(!is.na(Age_f), !is.na(paper_arm))

age_cols <- science_aaas_age_palette(age_levels)

fr <- if (file.exists(fr_path)) {
  read_tsv(fr_path, show_col_types = FALSE) |>
    mutate(paper_arm = factor(paper_arm, levels = arm_levels))
} else {
  tibble(
    paper_arm = factor(character(), levels = arm_levels),
    Friedman_p = numeric(),
    q_BH_Friedman_across_arms = numeric()
  )
}

if (nrow(fr) > 0L && "q_BH_Friedman_across_arms" %in% names(fr)) {
  fr <- fr |>
    transmute(
      paper_arm,
      fp = Friedman_p,
      fq = q_BH_Friedman_across_arms
    )
} else {
  fr <- tibble(
    paper_arm = factor(character(0L), levels = arm_levels),
    fp = numeric(0L),
    fq = numeric(0L)
  )
}

fmt_p <- function(p) {
  if (length(p) != 1L) {
    stop("fmt_p expects scalar")
  }
  if (is.na(p)) {
    return("\u2014")
  }
  if (p < 1e-3) {
    return(formatC(p, format = "e", digits = 1))
  }
  format(signif(p, 2), scientific = FALSE, trim = TRUE)
}

fr_lab <- fr |>
  rowwise() |>
  mutate(
    sig = !is.na(fq) && fq < 0.05,
    friedman_lab = if (is.na(fp) && is.na(fq)) {
      "Friedman (4 ages)\n\u2014"
    } else {
      paste0(
        sprintf(
          "Friedman (4 ages)\np = %s\nq_adj = %s",
          fmt_p(fp),
          fmt_p(fq)
        ),
        if (sig) "\n*" else ""
      )
    }
  ) |>
  ungroup() |>
  select(-sig)

y_ann <- df |>
  group_by(paper_arm) |>
  summarise(y_top = max(arg_frac, na.rm = TRUE), .groups = "drop") |>
  left_join(fr_lab, by = "paper_arm") |>
  mutate(
    friedman_lab = ifelse(is.na(friedman_lab), "Friedman (4 ages)\n\u2014", friedman_lab),
    y_label = y_top * 1.14,
    x_label = 2.5
  )

p <- ggplot(df, aes(x = Age_f, y = arg_frac, fill = Age_f)) +
  facet_wrap(~paper_arm, nrow = 2L, scales = "free_y") +
  geom_boxplot(
    width = 0.62,
    alpha = 0.45,
    colour = "#333333",
    linewidth = 0.4,
    outlier.shape = NA
  ) +
  geom_point(
    aes(color = Age_f),
    position = position_jitter(width = 0.12, height = 0, seed = 1L),
    size = 1.4,
    alpha = 0.82
  ) +
  geom_text(
    data = y_ann,
    aes(x = x_label, y = y_label, label = friedman_lab),
    inherit.aes = FALSE,
    size = 2.65,
    lineheight = 0.95,
    colour = "#1A1A1A",
    fontface = "plain"
  ) +
  scale_fill_manual(values = age_cols, name = "cohort_Age (d)", drop = FALSE) +
  scale_color_manual(values = age_cols, name = "cohort_Age (d)", drop = FALSE) +
  scale_x_discrete(name = "cohort_Age (d)") +
  scale_y_continuous(
    name = "ARG-carrying phage TPM / total viral TPM",
    expand = expansion(mult = c(0.04, 0.18))
  ) +
  labs(
    title = "PRJNA1010706: ARG phage fraction by treatment arm",
    subtitle = paste0(
      "Colour = sampling day (d). Friedman = pig-complete blocks across 4 ages; ",
      "q_adj = BH-FDR across arms. * : q_adj < 0.05."
    )
  ) +
  theme_science_aaas(11) +
  theme(
    axis.text.x = element_text(size = 9.5),
    legend.margin = margin(b = 0, t = 0)
  ) +
  guides(color = "none", fill = guide_legend(nrow = 1L, byrow = TRUE))

out_pdf <- file.path(fig_dir, "Fig_arg_frac_boxplot_lines_by_arm_facet.pdf")
out_png <- file.path(fig_dir, "Fig_arg_frac_boxplot_lines_by_arm_facet.png")

ggsave(
  out_pdf,
  p,
  width = 260 / 25.4,
  height = 150 / 25.4,
  device = "pdf"
)
ggsave(
  out_png,
  p,
  width = 260 / 25.4,
  height = 150 / 25.4,
  dpi = 320,
  bg = "white"
)

message("Wrote:\n  ", out_pdf, "\n  ", out_png)
