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
# FigB: heatmap of VIBRANT metabolism enrichment by host genus (log2FC) with FDR marks

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- args_all[grepl("^--file=", args_all)]
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else NA_character_
root <- file.path(PHAGE_ROOT, "Functional") else {
  normalizePath(".", winslash = "/")
}

in_tsv <- file.path(root, "summary_tables/vibrant_metabolism_host_enrichment/enrichment_metabolism_by_host_genus.tsv")
out_dir <- file.path(root, "figures/vibrant_metabolism_host_enrichment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

pal_hex <- c("#8FB4DC", "#FFDD8E", "#70CDBE", "#AC99D2", "#7AC3DF", "#F5AA61", "#EB756C")

stopifnot(file.exists(in_tsv))
df <- read_tsv(in_tsv, show_col_types = FALSE)
need <- c("Metabolism", "host_genus", "log2FC", "q_fdr")
stopifnot(all(need %in% names(df)))

df <- df %>%
  mutate(
    Metabolism = as.character(Metabolism),
    host_genus = as.character(host_genus),
    log2FC = as.numeric(log2FC),
    q_fdr = as.numeric(q_fdr),
    sig = ifelse(!is.na(q_fdr) & q_fdr < 0.05, "*", "")
  )

# choose a stable order: genera by total observed vOTUs across metabolism
if ("n_vOTUs_obs" %in% names(df)) {
  genus_order <- df %>%
    group_by(host_genus) %>%
    summarise(total_obs = sum(as.numeric(n_vOTUs_obs), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total_obs), host_genus) %>%
    pull(host_genus)
} else {
  genus_order <- sort(unique(df$host_genus))
}
metab_order <- sort(unique(df$Metabolism))

# complete full grid (40 genera x metabolism categories expected)
df2 <- df %>%
  select(Metabolism, host_genus, log2FC, q_fdr, sig) %>%
  complete(Metabolism = metab_order, host_genus = genus_order, fill = list(log2FC = 0, q_fdr = NA_real_, sig = ""))

# ggplot draws first factor level at bottom; reverse to make order top-to-bottom
df2$host_genus <- factor(df2$host_genus, levels = rev(genus_order))
df2$Metabolism <- factor(df2$Metabolism, levels = metab_order)

# diverging cool palette: green -> white -> purple
col_low <- pal_hex[3]
col_mid <- "white"
col_high <- pal_hex[4]

vmax <- max(abs(df2$log2FC), na.rm = TRUE)
vmax <- max(vmax, 0.5)

p <- ggplot(df2, aes(x = Metabolism, y = host_genus, fill = log2FC)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = sig), size = 5.2, color = "black") +
  scale_fill_gradient2(
    low = col_low,
    mid = col_mid,
    high = col_high,
    midpoint = 0,
    limits = c(-vmax, vmax),
    name = expression(log[2]*" fold-change")
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 13, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 14),
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.6, color = "gray25"),
    axis.ticks = element_line(linewidth = 0.5, color = "gray25"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.height = unit(18, "pt"),
    legend.key.width = unit(12, "pt"),
    legend.position = "top",
    legend.direction = "horizontal",
    plot.margin = margin(10, 10, 10, 10)
  )

# canvas size (scale with number of genera)
ng <- length(genus_order)
nm <- length(metab_order)
hh <- max(7.5, 0.22 * ng + 2.0)
ww <- max(8.5, 0.55 * nm + 4.0)

out_png <- file.path(out_dir, "FigB_heatmap_metabolism_by_host_with_FDR.png")
out_pdf <- file.path(out_dir, "FigB_heatmap_metabolism_by_host_with_FDR.pdf")
ggsave(out_png, p, width = ww, height = hh, dpi = 320)
ggsave(out_pdf, p, width = ww, height = hh, device = cairo_pdf)

## Top10 version (host genera)
top_n <- 10
genus_top <- head(genus_order, top_n)
# for the final figure, force Prevotella to the front (if present)
genus_top_final <- genus_top
if ("Prevotella" %in% genus_top_final) {
  genus_top_final <- c("Prevotella", genus_top_final[genus_top_final != "Prevotella"])
}

# Drop three low-informative metabolism columns in the Top10 figure (as requested)
drop_metab <- c(
  "Amino acid metabolism",
  "Folding, sorting and degradation",
  "Xenobiotics biodegradation and metabolism"
)

df2_top <- df2 %>%
  filter(as.character(host_genus) %in% genus_top) %>%
  filter(!(as.character(Metabolism) %in% drop_metab))
# keep top10 heatmap (non-final) as top-to-bottom
df2_top$host_genus <- factor(as.character(df2_top$host_genus), levels = rev(genus_top))
metab_top <- metab_order[!(metab_order %in% drop_metab)]
df2_top$Metabolism <- factor(as.character(df2_top$Metabolism), levels = metab_top)

ng2 <- length(genus_top)
hh2 <- max(6.0, 0.32 * ng2 + 2.0)
nm2 <- length(metab_top)
ww2 <- max(7.6, 0.55 * nm2 + 3.6)

p_top <- ggplot(df2_top, aes(x = Metabolism, y = host_genus, fill = log2FC)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = sig), size = 5.2, color = "black") +
  scale_fill_gradient2(
    low = col_low,
    mid = col_mid,
    high = col_high,
    midpoint = 0,
    limits = c(-vmax, vmax),
    name = expression(log[2]*" fold-change")
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 13, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 15),
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.6, color = "gray25"),
    axis.ticks = element_line(linewidth = 0.5, color = "gray25"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.height = unit(18, "pt"),
    legend.key.width = unit(12, "pt"),
    legend.position = "top",
    legend.direction = "horizontal",
    plot.margin = margin(10, 10, 10, 10)
  )

out_png2 <- file.path(out_dir, "FigB_heatmap_metabolism_by_host_with_FDR_top10.png")
out_pdf2 <- file.path(out_dir, "FigB_heatmap_metabolism_by_host_with_FDR_top10.pdf")
ggsave(out_png2, p_top, width = ww2, height = hh2, dpi = 320)
ggsave(out_pdf2, p_top, width = ww2, height = hh2, device = cairo_pdf)

## Top10 final version: swap x/y axes (x = host genus, y = Metabolism)
ww3 <- max(7.8, 0.50 * ng2 + 3.2)
hh3 <- max(6.2, 0.38 * nm2 + 2.2)

p_top_final_df <- df2_top
# for the final plot (x axis), enforce swapped left-to-right order
p_top_final_df$host_genus <- factor(as.character(p_top_final_df$host_genus), levels = genus_top_final)

p_top_final <- ggplot(p_top_final_df, aes(x = host_genus, y = Metabolism, fill = log2FC)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = sig), size = 5.2, color = "black") +
  scale_fill_gradient2(
    low = col_low,
    mid = col_mid,
    high = col_high,
    midpoint = 0,
    limits = c(-vmax, vmax),
    name = expression(log[2]*" fold-change")
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(size = 14, angle = 35, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 14),
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.6, color = "gray25"),
    axis.ticks = element_line(linewidth = 0.5, color = "gray25"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.height = unit(18, "pt"),
    legend.key.width = unit(12, "pt"),
    legend.position = "top",
    legend.direction = "horizontal",
    plot.margin = margin(10, 10, 10, 10)
  )

out_png3 <- file.path(out_dir, "FigB_heatmap_metabolism_by_host_with_FDR_top10_final.png")
out_pdf3 <- file.path(out_dir, "FigB_heatmap_metabolism_by_host_with_FDR_top10_final.pdf")
ggsave(out_png3, p_top_final, width = ww3, height = hh3, dpi = 320)
ggsave(out_pdf3, p_top_final, width = ww3, height = hh3, device = cairo_pdf)

message("Wrote:\n  ", out_png, "\n  ", out_pdf, "\n  ", out_png2, "\n  ", out_pdf2, "\n  ", out_png3, "\n  ", out_pdf3)

