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
# Redraw Fig_ge05_votu_host_genus_by_gut_stacked (R / ggplot2)
# Read meta_qc/ge05_votu_host_genus_tpm_by_gut.tsv (from analyze_ge05_votu_host_genus_tpm_by_gut.py).
# Colors: matplotlib tab20, order matches wide-table mean_pct_* columns and stack layers.
# Legend top, 3 rows; larger font; narrow geom_col width for bar spacing.
#
# Usage (project root):
#   Rscript meta_qc/plot_ge05_votu_host_genus_by_gut_stacked.R
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
if (!file.exists(file.path(root, "meta_qc", "ge05_votu_host_genus_tpm_by_gut.tsv"))) {
  root <- getwd()
}

OUT_DIR <- file.path(root, "meta_qc")
TSV_WIDE <- file.path(OUT_DIR, "ge05_votu_host_genus_tpm_by_gut.tsv")

GUT_ORDER <- c("Jejunum", "Ileum", "cecum", "Colon", "feces")

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

if (!file.exists(TSV_WIDE)) {
  stop("No ", TSV_WIDE, "；Please run: python analyze_ge05_votu_host_genus_tpm_by_gut.py")
}

wide <- fread(TSV_WIDE, sep = "\t")
ord <- match(GUT_ORDER, wide$Gut_harmonized)
wide <- wide[ord[!is.na(ord)]]

meas <- grep("^mean_pct_", names(wide), value = TRUE)
if (length(meas) < 1L) stop("table no mean_pct_* column")

lvl <- sub("^mean_pct_", "", meas)
plot_long <- melt(
  as.data.table(wide),
  id.vars = c("Gut_harmonized", "n_samples"),
  measure.vars = meas,
  variable.name = "var",
  value.name = "pct"
)
plot_long[, genus := factor(sub("^mean_pct_", "", as.character(var)), levels = lvl)]
plot_long <- plot_long[is.finite(pct)]
gut_lv <- gut_xlabel(wide$Gut_harmonized)
plot_long[, Gut := factor(gut_xlabel(Gut_harmonized), levels = gut_lv)]

# Genus → color fixed to match reference / Python tab20
GENUS_COLORS <- c(
  Unassigned = "#1f77b4",
  Prevotella = "#aec7e8",
  Lactobacillus = "#ff7f0e",
  Limosilactobacillus = "#ffbb78",
  Alloprevotella = "#2ca02c",
  Escherichia = "#98df8a",
  Gemmiger = "#d62728",
  Vescimonas = "#ff9896",
  Sodaliphilus = "#9467bd",
  Blautia_A = "#c5b0d5",
  Cryptobacteroides = "#8c564b",
  Faecousia = "#c49c94",
  Other = "#e377c2"
)
tab20_extra <- c(
  "#f7b6d2", "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5"
)
pal <- unname(GENUS_COLORS[lvl])
if (any(is.na(pal))) {
  miss <- which(is.na(pal))
  pal[miss] <- tab20_extra[(seq_along(miss) - 1L) %% length(tab20_extra) + 1L]
}
names(pal) <- lvl

LEGEND_TITLE <- "Host genus (prev>0.5 vOTUs)"

p <- ggplot(plot_long, aes(x = Gut, y = pct, fill = genus)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.28) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_fill_manual(
    values = pal,
    breaks = lvl,
    name = LEGEND_TITLE,
    guide = guide_legend(nrow = 3L, byrow = TRUE)
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.06, 0.06))) +
  labs(x = NULL, y = "Mean relative abundance (%)") +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    axis.text = element_text(size = 18, colour = "black"),
    axis.title.y = element_text(size = 20, margin = margin(r = 10)),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 16, face = "plain"),
    legend.key.height = unit(0.85, "lines"),
    legend.key.width = unit(1.1, "lines"),
    legend.spacing.x = unit(0.35, "lines"),
    legend.margin = margin(b = 4, t = 2),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(28, 8, 10, 6, "pt")
  )

ggsave(file.path(OUT_DIR, "Fig_ge05_votu_host_genus_by_gut_stacked.pdf"), p, width = 10.5, height = 7.2, useDingbats = FALSE)
ggsave(file.path(OUT_DIR, "Fig_ge05_votu_host_genus_by_gut_stacked.png"), p, width = 10.5, height = 7.2, dpi = 300)