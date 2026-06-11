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


library(ggplot2)
library(data.table)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- args_all[grepl("^--file=", args_all)]
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else NA_character_
root <- file.path(PHAGE_ROOT, "Functional") else {
  normalizePath(".", winslash = "/")
}

out_dir <- file.path(root, "figures/arg_country_threshold_compare")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

strict_in <- file.path(root, "summary_tables/arg_phage_carriers_sarg_only/sample_ARG_load.tsv")
nonstrict_in <- file.path(root, "summary_tables/arg_phage_carriers_sarg_nonstrict/sample_ARG_load.tsv")
cmp_in <- file.path(root, "summary_tables/arg_country_compare/country_ARG_load_compare_strict_vs_nonstrict.tsv")

df_s <- fread(strict_in)
df_n <- fread(nonstrict_in)

stopifnot(all(c("sample_id","Country","ARG_copyTPM_sum") %in% names(df_s)))
stopifnot(all(c("sample_id","Country","ARG_copyTPM_sum") %in% names(df_n)))

df_s[, threshold := "Strict (pident≥80, qcov≥0.85)"]
df_n[, threshold := "Non-strict (no pident/qcov filter)"]

df <- rbindlist(list(df_s, df_n), use.names = TRUE, fill = TRUE)
df[, ARG_copyTPM_sum := as.numeric(ARG_copyTPM_sum)]
df[, log10_load := log10(ARG_copyTPM_sum + 1e-6)]

# Country order: by strict median (descending) to keep interpretation stable
med_strict <- df[threshold == "Strict (pident≥80, qcov≥0.85)",
                 .(med = median(log10_load, na.rm = TRUE), n = .N), by = Country][order(-med)]
country_levels <- med_strict$Country
df[, Country := factor(Country, levels = rev(country_levels))]

# annotate sample size per country (based on strict table)
ann <- med_strict[, .(Country, n)]
ann[, Country := factor(Country, levels = country_levels)]

## Figure 1a: strict only
df_strict <- df[threshold == "Strict (pident≥80, qcov≥0.85)"]
p1a <- ggplot(df_strict, aes(x = Country, y = log10_load)) +
  geom_boxplot(
    width = 0.72,
    outlier.size = 0.25,
    outlier.alpha = 0.35,
    linewidth = 0.35,
    fill = "#70CDBE",
    color = "grey20"
  ) +
  geom_jitter(width = 0.12, size = 0.6, alpha = 0.10, color = "black") +
  coord_flip() +
  labs(
    x = NULL,
    y = expression(log[10]*"(ARG copy-TPM sum + 1e-6)"),
    title = "Country-level ARG burden (SARG strict)"
  ) +
  theme_classic(base_size = 22) +
  theme(
    text = element_text(color = "black"),
    axis.text.y = element_text(size = 16.5, color = "black"),
    axis.text.x = element_text(size = 16.5, color = "black"),
    axis.title.x = element_text(size = 18.5, color = "black"),
    plot.title = element_text(face = "bold", size = 22, color = "black"),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.35),
    panel.grid.major.y = element_blank()
  )

ggsave(file.path(out_dir, "Fig1a_box_country_ARGload_SARG_strict.pdf"),
       p1a, width = 5.8, height = 9.0, device = cairo_pdf)
ggsave(file.path(out_dir, "Fig1a_box_country_ARGload_SARG_strict.png"),
       p1a, width = 5.8, height = 9.0, dpi = 450)

## Figure 1a (mean): strict only, country mean ± 95% CI on raw scale (ARG_copyTPM_sum)
mean_tab <- df_strict[, .(
  n = .N,
  mean_load = mean(ARG_copyTPM_sum, na.rm = TRUE),
  sd_load = sd(ARG_copyTPM_sum, na.rm = TRUE)
), by = Country]
mean_tab[, se_load := sd_load / sqrt(pmax(1, n))]
mean_tab[, ci_lo := mean_load - 1.96 * se_load]
mean_tab[, ci_hi := mean_load + 1.96 * se_load]
mean_tab[ci_lo < 0, ci_lo := 0]

# order by mean desc
mean_tab <- mean_tab[order(-mean_load)]
mean_tab[, Country_lab := paste0(Country, " (n=", n, ")")]
mean_tab[, Country_lab := factor(Country_lab, levels = rev(mean_tab$Country_lab))]

p1a_mean <- ggplot(mean_tab, aes(x = Country_lab, y = mean_load)) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.25, linewidth = 0.7, color = "#C7B9E8") +
  geom_point(size = 3.5, color = "#AC99D2") +
  coord_flip() +
  labs(
    x = NULL,
    y = "Mean ARG copy-TPM sum (SARG strict)",
    title = "Country-level mean ARG burden (SARG strict)"
  ) +
  theme_classic(base_size = 22) +
  theme(
    text = element_text(color = "black"),
    axis.text.y = element_text(size = 15.5, color = "black"),
    axis.text.x = element_text(size = 15.5, color = "black"),
    axis.title.x = element_text(size = 18.5, color = "black"),
    plot.title = element_text(face = "bold", size = 22, color = "black"),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.35),
    panel.grid.major.y = element_blank()
  )

fwrite(mean_tab[, .(Country, n, mean_load, sd_load, se_load, ci_lo, ci_hi)],
       file.path(out_dir, "Fig1a_country_mean_ARGload_SARG_strict.tsv"),
       sep = "\t")
ggsave(file.path(out_dir, "Fig1a_mean_country_ARGload_SARG_strict.pdf"),
       p1a_mean, width = 7.2, height = 9.2, device = cairo_pdf)
ggsave(file.path(out_dir, "Fig1a_mean_country_ARGload_SARG_strict.png"),
       p1a_mean, width = 7.2, height = 9.2, dpi = 450)

## Figure 1b: non-strict only
df_nonstrict <- df[threshold == "Non-strict (no pident/qcov filter)"]
p1b <- ggplot(df_nonstrict, aes(x = Country, y = log10_load)) +
  geom_boxplot(
    width = 0.72,
    outlier.size = 0.25,
    outlier.alpha = 0.35,
    linewidth = 0.35,
    fill = "#70CDBE",
    color = "grey20"
  ) +
  geom_jitter(width = 0.12, size = 0.6, alpha = 0.10, color = "black") +
  coord_flip() +
  labs(
    x = NULL,
    y = expression(log[10]*"(ARG copy-TPM sum + 1e-6)"),
    title = "Country-level ARG burden (SARG non-strict)"
  ) +
  theme_classic(base_size = 22) +
  theme(
    text = element_text(color = "black"),
    axis.text.y = element_text(size = 16.5, color = "black"),
    axis.text.x = element_text(size = 16.5, color = "black"),
    axis.title.x = element_text(size = 18.5, color = "black"),
    plot.title = element_text(face = "bold", size = 22, color = "black"),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.35),
    panel.grid.major.y = element_blank()
  )

ggsave(file.path(out_dir, "Fig1b_box_country_ARGload_SARG_nonstrict.pdf"),
       p1b, width = 5.8, height = 9.0, device = cairo_pdf)
ggsave(file.path(out_dir, "Fig1b_box_country_ARGload_SARG_nonstrict.png"),
       p1b, width = 5.8, height = 9.0, dpi = 450)

## Figure 1b (mean): non-strict, country mean ± 95% CI on raw scale (ARG_copyTPM_sum)
mean_tab_n <- df_nonstrict[, .(
  n = .N,
  mean_load = mean(ARG_copyTPM_sum, na.rm = TRUE),
  sd_load = sd(ARG_copyTPM_sum, na.rm = TRUE)
), by = Country]
mean_tab_n[, se_load := sd_load / sqrt(pmax(1, n))]
mean_tab_n[, ci_lo := mean_load - 1.96 * se_load]
mean_tab_n[, ci_hi := mean_load + 1.96 * se_load]
mean_tab_n[ci_lo < 0, ci_lo := 0]

# order by mean desc
mean_tab_n <- mean_tab_n[order(-mean_load)]
mean_tab_n[, Country_lab := paste0(Country, " (n=", n, ")")]
mean_tab_n[, Country_lab := factor(Country_lab, levels = rev(mean_tab_n$Country_lab))]

p1b_mean <- ggplot(mean_tab_n, aes(x = Country_lab, y = mean_load)) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.25, linewidth = 0.7, color = "#C7B9E8") +
  geom_point(size = 3.5, color = "#AC99D2") +
  coord_flip() +
  labs(
    x = NULL,
    y = "Mean ARG copy-TPM sum (SARG non-strict)",
    title = "Country-level mean ARG burden (SARG non-strict)"
  ) +
  theme_classic(base_size = 22) +
  theme(
    text = element_text(color = "black"),
    axis.text.y = element_text(size = 15.5, color = "black"),
    axis.text.x = element_text(size = 15.5, color = "black"),
    axis.title.x = element_text(size = 18.5, color = "black"),
    plot.title = element_text(face = "bold", size = 22, color = "black"),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.35),
    panel.grid.major.y = element_blank()
  )

fwrite(mean_tab_n[, .(Country, n, mean_load, sd_load, se_load, ci_lo, ci_hi)],
       file.path(out_dir, "Fig1b_country_mean_ARGload_SARG_nonstrict.tsv"),
       sep = "\t")
ggsave(file.path(out_dir, "Fig1b_mean_country_ARGload_SARG_nonstrict.pdf"),
       p1b_mean, width = 7.2, height = 9.2, device = cairo_pdf)
ggsave(file.path(out_dir, "Fig1b_mean_country_ARGload_SARG_nonstrict.png"),
       p1b_mean, width = 7.2, height = 9.2, dpi = 450)

## Figure 2: per-country effect size (median fold-change)
cmp <- fread(cmp_in)
stopifnot(all(c("Country","ARG_copyTPM_sum_median_strict","ARG_copyTPM_sum_median_nonstrict","n_samples_strict") %in% names(cmp)))
cmp[, median_strict := as.numeric(ARG_copyTPM_sum_median_strict)]
cmp[, median_nonstrict := as.numeric(ARG_copyTPM_sum_median_nonstrict)]
cmp[, n := as.integer(n_samples_strict)]
cmp[, log2_fc := log2((median_nonstrict + 1e-12) / (median_strict + 1e-12))]

# order by log2 fold-change
cmp <- cmp[order(log2_fc)]
cmp[, Country := factor(Country, levels = cmp$Country)]

p2 <- ggplot(cmp, aes(x = Country, y = log2_fc)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "grey45") +
  geom_segment(aes(xend = Country, y = 0, yend = log2_fc),
               linewidth = 1.0, color = "grey70") +
  geom_point(aes(size = n, color = log2_fc),
             alpha = 0.95) +
  coord_flip() +
  scale_color_gradient2(low = "#2C7FB8", mid = "grey80", high = "#F03B20", midpoint = 0) +
  scale_size_continuous(range = c(2.0, 8.0)) +
  labs(
    x = NULL,
    y = expression(log[2]*"(median non-strict / median strict)"),
    color = "log2 FC",
    size = "N samples",
    title = "Effect of threshold relaxation on country-level ARG burden"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 15)
  )

ggsave(file.path(out_dir, "Fig2_lollipop_country_log2FC_nonstrict_vs_strict.pdf"),
       p2, width = 10.8, height = 9.0, device = cairo_pdf)
ggsave(file.path(out_dir, "Fig2_lollipop_country_log2FC_nonstrict_vs_strict.png"),
       p2, width = 10.8, height = 9.0, dpi = 360)

