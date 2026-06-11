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
# Provirus maps: linear + circular (gggenes, VIBRANT legend, kb scale).

suppressPackageStartupMessages({
  library(ggplot2)
  library(gggenes)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg)) {
  root <- dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/"))
} else {
  root <- normalizePath(getwd(), winslash = "/")
}

source(file.path(root, "R", "science_aaas_figures.R"), local = TRUE)
source(file.path(root, "R", "provirus_func_palette.R"), local = TRUE)

gene_tsv <- file.path(
  root, "results", "phage_host_species_phenotype", "figures",
  "SRR15732359_provirus_147_40920_gggenes_genes.tsv"
)
out_dir <- dirname(gene_tsv)

genes <- read.delim(gene_tsv, stringsAsFactors = FALSE)
genes$strand <- ifelse(genes$strand >= 0, 1L, -1L)
if (!"plot_label" %in% names(genes)) genes$plot_label <- ""
genes$plot_label[is.na(genes$plot_label)] <- ""
genes$func_category <- factor(genes$func_category, levels = provirus_func_levels())

genome_len <- max(genes$end, na.rm = TRUE)
pal <- provirus_func_palette()
lab_func <- provirus_func_labels()
func_levels <- provirus_func_levels()

out_pdf <- file.path(out_dir, "SRR15732359_provirus_linear.pdf")
out_png <- file.path(out_dir, "SRR15732359_provirus_linear.png")
out_pdf_clean <- file.path(out_dir, "SRR15732359_provirus_linear_nolabel.pdf")
out_png_clean <- file.path(out_dir, "SRR15732359_provirus_linear_nolabel.png")
out_circ_pdf <- file.path(out_dir, "SRR15732359_provirus_circular.pdf")
out_circ_png <- file.path(out_dir, "SRR15732359_provirus_circular.png")
out_circ_pdf_clean <- file.path(out_dir, "SRR15732359_provirus_circular_nolabel.pdf")
out_circ_png_clean <- file.path(out_dir, "SRR15732359_provirus_circular_nolabel.png")

genes$molecule <- "provirus"
genes$forward <- genes$strand >= 0L

present <- intersect(func_levels, unique(as.character(genes$func_category)))

short_label <- function(txt) {
  if (!nzchar(txt)) return("")
  x <- gsub('["\']', "", txt)
  x <- sub("^REFSEQ\\s+", "", x, ignore.case = TRUE)
  x <- sub("^sp\\|[^|]+\\|", "", x)
  x <- sub("^K[0-9]{5};\\s*", "", x)
  x <- sub(";.*$", "", x)
  x <- sub("\\s+protein$", "", x, ignore.case = TRUE)
  if (nchar(x) > 18) x <- paste0(substr(x, 1, 15), "...")
  x
}

pick_labels <- function(g, n = 10L) {
  if (n <= 0L) return(g[0, , drop = FALSE])
  pri <- c(
    packaging = 1L, tail = 2L, lysogeny = 3L, lysis = 4L,
    stability = 5L, metabolism = 6L, dna_replication = 7L,
    regulation = 8L
  )
  g$lab <- vapply(seq_len(nrow(g)), function(i) short_label(g$plot_label[i]), "")
  c <- g[nzchar(g$lab) & g$func_category != "hypothetical", , drop = FALSE]
  if (nrow(c) == 0) return(c)
  c$pri <- pri[as.character(c$func_category)]
  c$pri[is.na(c$pri)] <- 50L
  c$glen <- c$end - c$start + 1L
  c <- c[order(c$pri, -c$glen), , drop = FALSE]
  out <- c[FALSE, , drop = FALSE]
  used <- character()
  for (i in seq_len(nrow(c))) {
    if (c$lab[i] %in% used) next
    used <- c(used, c$lab[i])
    out <- rbind(out, c[i, , drop = FALSE])
    if (nrow(out) >= n) break
  }
  out
}

caption <- sprintf(
  "SRR15732359_viral_contig_2826|provirus_147_40920 (%s bp, 54 ORFs)",
  format(genome_len, big.mark = ",")
)

x_breaks <- c(seq(0L, 40000L, by = 5000L), genome_len)

format_kb <- function(x) {
  ifelse(
    abs(x - genome_len) < 300,
    sprintf("%.1f", x / 1000),
    as.character(as.integer(x / 1000))
  )
}

base_plot <- function() {
  ggplot(genes, aes(
    xmin = start, xmax = end, y = molecule,
    fill = func_category, forward = forward
  )) +
    geom_gene_arrow(
      arrowhead_height = grid::unit(2.4, "mm"),
      arrowhead_width = grid::unit(2, "mm"),
      arrow_body_height = grid::unit(2.6, "mm"),
      colour = "black",
      linewidth = 0.28
    ) +
    scale_fill_manual(
      name = "Function (VIBRANT)",
      values = pal[present],
      labels = lab_func[present],
      breaks = present,
      drop = FALSE
    ) +
    guides(
      fill = guide_legend(
        ncol = 2,
        byrow = TRUE,
        keywidth = grid::unit(0.45, "cm"),
        keyheight = grid::unit(0.32, "cm"),
        title.position = "top",
        title.hjust = 0.5
      )
    )
}

add_gene_labels <- function(p, gl) {
  if (nrow(gl) == 0) return(p)
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p + ggrepel::geom_text_repel(
      data = gl,
      aes(x = x, y = molecule, label = lab),
      inherit.aes = FALSE,
      size = 2.7, fontface = "italic", colour = "black",
      min.segment.length = 0.1, segment.size = 0.22, segment.colour = "grey35",
      box.padding = 0.25, point.padding = 0.15, max.overlaps = 25,
      force = 3, force_pull = 0.5, direction = "both"
    )
  } else {
    p + geom_text(
      data = gl,
      aes(x = x, y = molecule, label = lab),
      inherit.aes = FALSE,
      size = 2.5, fontface = "italic", vjust = -0.8
    )
  }
}

plot_linear <- function(show_gene_labels = TRUE) {
  gl <- pick_labels(genes, if (show_gene_labels) 10L else 0L)
  if (nrow(gl) > 0) {
    gl$x <- (gl$start + gl$end) / 2
    gl$lab <- vapply(seq_len(nrow(gl)), function(i) short_label(gl$plot_label[i]), "")
  }

  p <- base_plot() +
    scale_x_continuous(
      name = "Genomic position (kb)",
      limits = c(-500, genome_len + 500),
      breaks = x_breaks,
      labels = format_kb,
      minor_breaks = seq(0L, genome_len, by = 1000L),
      expand = c(0, 0)
    ) +
    labs(
      x = "Genomic position (kb)", y = NULL, title = NULL,
      caption = paste0(
        caption,
        "  |  Arrow: transcription 5'-3'; right = plus strand, left = minus strand"
      )
    ) +
    theme_genes() +
    theme(
      plot.margin = margin(8, 12, 8, 12),
      legend.position = "bottom",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 7.5),
      legend.box.spacing = grid::unit(0.3, "cm"),
      legend.margin = margin(t = 2, b = 2),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 7.5, colour = "grey15"),
      axis.ticks.x = element_line(colour = "grey30", linewidth = 0.35),
      axis.title.x = element_text(size = 8.5, colour = "grey15", margin = margin(t = 6)),
      axis.line.x = element_line(colour = "grey35", linewidth = 0.4),
      panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.25),
      panel.grid.minor.x = element_line(colour = "grey95", linewidth = 0.15),
      axis.ticks.length.x = grid::unit(2.5, "pt"),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.caption = element_text(size = 7.5, colour = "grey25", hjust = 0.5, margin = margin(t = 6))
    ) +
    coord_cartesian(clip = "off") +
    annotate(
      "segment", x = 0, xend = 0, y = 0.92, yend = 1.08,
      linewidth = 0.55, colour = "black"
    ) +
    annotate(
      "segment", x = genome_len, xend = genome_len, y = 0.92, yend = 1.08,
      linewidth = 0.55, colour = "black"
    )

  add_gene_labels(p, gl)
}

plot_circular <- function(show_gene_labels = TRUE) {
  gl <- pick_labels(genes, if (show_gene_labels) 8L else 0L)
  if (nrow(gl) > 0) {
    gl$x <- (gl$start + gl$end) / 2
    gl$lab <- vapply(seq_len(nrow(gl)), function(i) short_label(gl$plot_label[i]), "")
  }

  p <- base_plot() +
    coord_polar(theta = "x", start = 0.25, direction = 1, clip = "off") +
    scale_x_continuous(
      limits = c(0, genome_len),
      breaks = x_breaks,
      labels = format_kb,
      minor_breaks = seq(0L, genome_len, by = 2000L),
      expand = c(0, 0)
    ) +
    labs(
      title = "SRR15732359 provirus",
      subtitle = paste0(
        caption,
        " | Circular | kb clockwise from 0 at 12 o'clock | ",
        "outward=plus, inward=minus"
      ),
      x = NULL,
      y = NULL
    ) +
    theme_genes() +
    theme(
      plot.margin = margin(2, 2, 6, 2),
      legend.position = "bottom",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 7.5),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 6.5, colour = "grey20"),
      axis.ticks.x = element_line(colour = "grey35", linewidth = 0.3),
      panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.2),
      panel.grid.minor.x = element_line(colour = "grey93", linewidth = 0.12),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 6.8, colour = "grey35", hjust = 0.5, lineheight = 0.95)
    )

  add_gene_labels(p, gl)
}

save_linear <- function(show_gene_labels, pdf_path, png_path) {
  h <- if (show_gene_labels) 3.9 else 3.5
  pdf(pdf_path, width = 14, height = h)
  print(plot_linear(show_gene_labels))
  dev.off()
  png(png_path, width = 14, height = h, units = "in", res = 300)
  print(plot_linear(show_gene_labels))
  dev.off()
}

save_circular <- function(show_gene_labels, pdf_path, png_path) {
  pdf(pdf_path, width = 9, height = 9)
  print(plot_circular(show_gene_labels))
  dev.off()
  png(png_path, width = 9, height = 9, units = "in", res = 300)
  print(plot_circular(show_gene_labels))
  dev.off()
}

save_linear(TRUE, out_pdf, out_png)
save_linear(FALSE, out_pdf_clean, out_png_clean)
save_circular(TRUE, out_circ_pdf, out_circ_png)
save_circular(FALSE, out_circ_pdf_clean, out_circ_png_clean)

message("Linear (labels)    -> ", out_pdf)
message("Linear (no labels) -> ", out_pdf_clean)
message("Circular (labels)  -> ", out_circ_pdf)
message("Circular (no lbl)  -> ", out_circ_pdf_clean)
