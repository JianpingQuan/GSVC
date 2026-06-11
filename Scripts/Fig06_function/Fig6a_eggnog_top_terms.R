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
# VIBRANT AMG distribution by KEGG Metabolism class — two separate figures (gene-level vs vOTU-level).
# Inputs: summary_tables/vibrant_metabolism_host_enrichment/overall_metabolism_*.tsv
# Outputs: figures/vibrant_metabolism_host_enrichment/Fig6b_* and FigS1_* (PNG + PDF)

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- args_all[grepl("^--file=", args_all)]
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else NA_character_
root <- file.path(PHAGE_ROOT, "Functional") else {
  normalizePath(".", winslash = "/")
}

in_genes <- file.path(root, "summary_tables/vibrant_metabolism_host_enrichment/overall_metabolism_gene_counts.tsv")
in_votu  <- file.path(root, "summary_tables/vibrant_metabolism_host_enrichment/overall_metabolism_votu_counts.tsv")
out_dir  <- file.path(root, "figures/vibrant_metabolism_host_enrichment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# Green → purple gradient (aligned with host-metabolism heatmap / project palette)
grad_low  <- "#BFE9C8"
grad_mid  <- "#9BC9D8"
grad_high <- "#6B4A9E"

plot_metabolism_bars <- function(df, value_lab, title, subtitle = NULL) {
  stopifnot(all(c("Metabolism", "n") %in% names(df)))

  df <- df %>%
    mutate(
      Metabolism = as.character(Metabolism),
      n = as.integer(n)
    ) %>%
    filter(!is.na(Metabolism), !is.na(n), n > 0L) %>%
    arrange(n)

  # Largest n at top: levels in ascending n order (first level = bottom in ggplot)
  df$Metabolism <- factor(df$Metabolism, levels = unique(df$Metabolism))

  ggplot(df, aes(x = n, y = Metabolism, fill = n)) +
    geom_col(width = 0.66, color = "white", linewidth = 0.22) +
    geom_text(
      aes(label = comma(n)),
      hjust = -0.08,
      size = 3.5,
      color = "gray15"
    ) +
    scale_fill_gradientn(
      colours = c(grad_low, grad_mid, grad_high),
      guide = "none"
    ) +
    scale_x_continuous(
      name = value_lab,
      expand = expansion(mult = c(0.02, 0.16)),
      labels = label_comma()
    ) +
    labs(title = title, subtitle = subtitle, y = NULL) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13.5),
      plot.subtitle = element_text(color = "gray35", size = 10),
      axis.title.y = element_blank(),
      axis.text.y = element_text(color = "gray10", size = 10.5),
      axis.text.x = element_text(color = "gray10", size = 10),
      axis.title.x = element_text(size = 11.5),
      axis.line = element_line(linewidth = 0.55, color = "gray20"),
      axis.ticks = element_line(linewidth = 0.5, color = "gray20"),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.4),
      panel.grid.major.y = element_blank()
    )
}

df_g <- read_tsv(in_genes, show_col_types = FALSE)
df_v <- read_tsv(in_votu, show_col_types = FALSE)

stopifnot(ncol(df_g) >= 2L, ncol(df_v) >= 2L)
df_g <- df_g %>% rename(Metabolism = 1, n = 2)
df_v <- df_v %>% rename(Metabolism = 1, n = 2)

p_gene <- plot_metabolism_bars(
  df_g,
  "Unique AMG genes (proteins)",
  "VIBRANT AMG distribution by metabolism (gene level)",
  subtitle = "Source: overall_metabolism_gene_counts.tsv"
)

p_votu <- plot_metabolism_bars(
  df_v,
  "vOTUs carrying \u2265 1 AMG in category",
  "VIBRANT AMG distribution by metabolism (vOTU level)",
  subtitle = "Source: overall_metabolism_votu_counts.tsv"
)

# Figure heights scale with number of categories
ng <- nrow(df_g)
nv <- nrow(df_v)
h_g <- max(4.8, 0.32 * ng + 1.2)
h_v <- max(4.8, 0.32 * nv + 1.2)
w_in <- 7.4

# Manuscript-oriented filenames: main panel vOTU = Fig 6b; supplement gene = Fig S1
out_g_png <- file.path(out_dir, "FigS1_vibrant_AMG_metabolism_genes.png")
out_g_pdf <- file.path(out_dir, "FigS1_vibrant_AMG_metabolism_genes.pdf")
out_v_png <- file.path(out_dir, "Fig6b_vibrant_AMG_metabolism_vOTU.png")
out_v_pdf <- file.path(out_dir, "Fig6b_vibrant_AMG_metabolism_vOTU.pdf")

ggsave(out_g_png, p_gene, width = w_in, height = h_g, dpi = 320)
ggsave(out_g_pdf, p_gene, width = w_in, height = h_g, device = cairo_pdf)
ggsave(out_v_png, p_votu, width = w_in, height = h_v, dpi = 320)
ggsave(out_v_pdf, p_votu, width = w_in, height = h_v, device = cairo_pdf)

message("Wrote:\n  ", out_g_png, "\n  ", out_g_pdf, "\n  ", out_v_png, "\n  ", out_v_pdf)
