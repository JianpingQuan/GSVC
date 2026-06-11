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
# Fig.3-style three-panel figure (from your TSV):
#   (a) Class-level donut (subset with completeness > threshold)
#   (b) Top: Family composition pie within named families only; bottom: Assigned_F vs Unclassified_F pie
#   (c) Class–Family chord diagram (circlize::chordDiagram)
#
# Main table must contain: seq_name, taxonomy (geNomad semicolon path, 6 segments after removing Viruses: … class … family)
# completeness optional:
#   - If main table has completeness column, use it directly;
#   - Otherwise merge via --completeness_tsv (default join on seq_name).
# If main table lacks completeness column but rows are already ">50% pre-filtered" (e.g. high/medium-quality representative set), add:
#   --assume_completeness_gt50
#   so caption and text are consistent (this script does not drop rows for this reason).
#
# Dependencies: readr, dplyr, tidyr, ggplot2, patchwork, scales, circlize, grDevices
#
# Rscript scripts/plot_fig3_taxonomy_composition_R.R \
#   --tsv path/high_medium_quality_viral_rep_seq_virus_summary.tsv \
#   --out_dir path/figures_nm
# Rscript ... --completeness_tsv path/checkv_quality_summary.tsv --min_completeness 50
# Rscript ... --top_n 10   # (a)(b top)(c) show Top N classes/families by count; rest merged into Other (default 10)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(circlize)
  library(grid)
})

OTHER_LAB <- "Other"

# Discrete colors similar to literature/NM-style legends (not bound to taxon names; assigned in legend order)
PALETTE_RANKED <- c(
  "#1B7837", "#312783", "#80CCE1", "#4DA699", "#D9C985",
  "#C96D7D", "#A64D9B", "#8C9433", "#6699CC", "#44AA99"
)
COL_UNASSIGNED <- "#B3B3B3"
COL_OTHER <- "#661400"

# levels_ordered: factor level order (e.g. decreasing frequency, Other last)
assign_fig3_fill_colors <- function(levels_ordered, mode = c("class", "family")) {
  mode <- match.arg(mode)
  out <- setNames(character(length(levels_ordered)), levels_ordered)
  j <- 1L
  np <- length(PALETTE_RANKED)
  for (nm in levels_ordered) {
    if (identical(nm, OTHER_LAB)) {
      out[[nm]] <- COL_OTHER
    } else if (mode == "class" && identical(nm, "Unclassified_C")) {
      out[[nm]] <- COL_UNASSIGNED
    } else if (mode == "family" && identical(nm, "Unclassified_F")) {
      out[[nm]] <- COL_UNASSIGNED
    } else {
      out[[nm]] <- PALETTE_RANKED[((j - 1L) %% np) + 1L]
      j <- j + 1L
    }
  }
  out
}

# Legend: ` 93.60% 4,960 Name` (percentage and count columns roughly aligned; with mono legend font)
legend_pct_count_name <- function(names_chr, pct, n) {
  n_fmt <- format(as.integer(n), big.mark = ",", trim = TRUE)
  nw <- max(nchar(n_fmt), 1L)
  n_pad <- sprintf(paste0("%", nw, "s"), n_fmt)
  pct_str <- sprintf("%6.2f%%", as.numeric(pct))
  paste(pct_str, n_pad, names_chr)
}

legend_labels_for_levels <- function(tab, level_vec, name_col) {
  nm <- as.character(tab[[name_col]])
  idx <- match(level_vec, nm)
  if (anyNA(idx)) {
    stop("legend match NA for ", name_col, call. = FALSE)
  }
  legend_pct_count_name(nm[idx], tab$pct[idx], tab$n[idx])
}

# Keep Top n by frequency; merge rest into Other (no Other if unique classes <= n)
collapse_top_other <- function(x, n, other_label = OTHER_LAB) {
  x <- as.character(x)
  tab <- sort(table(x), decreasing = TRUE)
  if (length(tab) <= n) {
    return(x)
  }
  keep <- names(tab)[seq_len(n)]
  ifelse(x %in% keep, x, other_label)
}

# Plot factor order: descending n, Other last
order_factor_other_last <- function(fct_chr, counts_named) {
  ord <- names(sort(counts_named, decreasing = TRUE))
  if (OTHER_LAB %in% ord) {
    ord <- c(setdiff(ord, OTHER_LAB), OTHER_LAB)
  }
  factor(fct_chr, levels = ord)
}

parse_class_family <- function(tax) {
  if (is.na(tax) || !nzchar(trimws(as.character(tax)))) {
    return(c(class = "Unclassified_C", family = "Unclassified_F"))
  }
  parts <- strsplit(trimws(as.character(tax)), ";", fixed = TRUE)[[1]]
  if (length(parts) && parts[[1]] == "Viruses") {
    parts <- parts[-1]
  }
  while (length(parts) < 6) {
    parts <- c(parts, "")
  }
  p <- trimws(parts[seq_len(6)])
  cls <- if (nzchar(p[[4]])) p[[4]] else "Unclassified_C"
  fam <- if (nzchar(p[[6]])) p[[6]] else "Unclassified_F"
  c(class = cls, family = fam)
}

get_arg <- function(argv, name, default) {
  i <- match(name, argv)
  if (is.na(i) || i >= length(argv)) {
    return(default)
  }
  argv[[i + 1]]
}

argv <- commandArgs(trailingOnly = TRUE)
tsv <- get_arg(argv, "--tsv", ""/high_medium_quality_viral_rep_seq_virus_summary.tsv")
out_dir <- get_arg(argv, "--out_dir", file.path(dirname(tsv), "..", "Functional", "figures_nm"))
out_base <- get_arg(argv, "--out_base", file.path(out_dir, "Fig3_taxonomy_composition_user_data"))
comp_tsv <- get_arg(argv, "--completeness_tsv", NA_character_)
comp_id <- get_arg(argv, "--completeness_id_col", "contig_id")
comp_col <- get_arg(argv, "--completeness_col", "completeness")
id_main <- get_arg(argv, "--id_col", "seq_name")
min_comp <- as.numeric(get_arg(argv, "--min_completeness", "50"))
top_n <- as.integer(get_arg(argv, "--top_n", "10"))
if (is.na(top_n) || top_n < 1L) {
  top_n <- 10L
}
assume_gt50 <- ("--assume_completeness_gt50" %in% argv) ||
  ("--assume-completeness-gt50" %in% argv)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

need_cols <- c(id_main, "taxonomy")
raw <- read_tsv(tsv, show_col_types = FALSE)
miss <- setdiff(need_cols, names(raw))
if (length(miss)) {
  stop("TSV less column: ", paste(miss, collapse = ", "), call. = FALSE)
}

comp_note <- character(0)
if (!is.na(comp_tsv) && nzchar(comp_tsv) && file.exists(comp_tsv)) {
  cv <- read_tsv(comp_tsv, show_col_types = FALSE)
  if (!comp_id %in% names(cv) || !comp_col %in% names(cv)) {
    stop("completeness less column: ", comp_id, " / ", comp_col, call. = FALSE)
  }
  cv <- cv %>%
    transmute(
      !!id_main := as.character(.data[[comp_id]]),
      .comp = suppressWarnings(as.numeric(.data[[comp_col]]))
    ) %>%
    filter(!is.na(.comp))
  raw <- raw %>%
    left_join(cv, by = id_main)
  raw <- raw %>% filter(.comp >= min_comp)
  comp_note <- sprintf("completeness >= %.1f%%（来自 %s）", min_comp, basename(comp_tsv))
} else if ("completeness" %in% names(raw)) {
  raw <- raw %>%
    mutate(.comp = suppressWarnings(as.numeric(completeness))) %>%
    filter(!is.na(.comp), .comp >= min_comp)
  comp_note <- sprintf("completeness >= %.1f%%（main table completeness column）", min_comp)
} else if (assume_gt50) {
  comp_note <- sprintf(
    paste0(
      "Main table no completeness column；this table have filted according to completeness > %.0f%%"
    ),
    min_comp, basename(tsv)
  )
} else {
  comp_note <- sprintf(
    "no completeness column and no --assume_completeness_gt50",
    basename(tsv)
  )
}

cf <- t(vapply(raw$taxonomy, parse_class_family, character(2)))
raw$class_lab <- cf[, 1]
raw$family_lab <- cf[, 2]
raw$assigned_f <- raw$family_lab != "Unclassified_F"

raw$class_plot <- collapse_top_other(raw$class_lab, top_n)
raw$family_plot <- collapse_top_other(raw$family_lab, top_n)

# ---- (a) Class donut（Top N + Other）----
tab_a <- raw %>%
  count(class_plot, name = "n") %>%
  mutate(pct = n / sum(n) * 100, frac = n / sum(n))
cn <- setNames(tab_a$n, tab_a$class_plot)
tab_a$class_plot <- order_factor_other_last(tab_a$class_plot, cn)

fills_a <- assign_fig3_fill_colors(levels(tab_a$class_plot), "class")
lv_a <- levels(tab_a$class_plot)
lab_a <- legend_labels_for_levels(tab_a, lv_a, "class_plot")

p_a <- ggplot(tab_a, aes(x = 2.6, y = frac, fill = class_plot)) +
  geom_col(width = 0.9, colour = "white", linewidth = 0.35) +
  coord_polar(theta = "y") +
  xlim(0.5, 3.2) +
  scale_fill_manual(
    values = fills_a,
    breaks = lv_a,
    labels = lab_a,
    name = NULL,
    drop = FALSE
  ) +
  theme_void(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 11, margin = margin(b = 4)),
    legend.position = "right",
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size = 7.5, family = "mono"),
    plot.background = element_rect(fill = "white", colour = NA)
  ) +
  labs(title = sprintf("(a) Class-level composition (top %d + Other)", top_n))

# ---- (b) pies ----
df_ass <- raw %>% filter(assigned_f)
if (nrow(df_ass) == 0L) {
  tab_f <- tibble(
    family_plot = factor("Unclassified_F", levels = "Unclassified_F"),
    n = nrow(raw),
    pct = 100,
    frac = 1
  )
} else {
  tab_f <- df_ass %>%
    count(family_plot, name = "n") %>%
    mutate(pct = n / sum(n) * 100, frac = n / sum(n))
  fn <- setNames(tab_f$n, tab_f$family_plot)
  tab_f$family_plot <- order_factor_other_last(tab_f$family_plot, fn)
}

pal_f <- assign_fig3_fill_colors(levels(tab_f$family_plot), "family")
lv_f <- levels(tab_f$family_plot)
lab_f <- legend_labels_for_levels(tab_f, lv_f, "family_plot")

p_b1 <- ggplot(tab_f, aes(x = "", y = frac, fill = family_plot)) +
  geom_col(width = 1, colour = "white", linewidth = 0.3) +
  coord_polar(theta = "y") +
  scale_fill_manual(
    values = pal_f,
    breaks = lv_f,
    labels = lab_f,
    name = NULL,
    drop = FALSE
  ) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 9.5, hjust = 0.5),
    legend.position = "right",
    legend.key.size = unit(0.28, "cm"),
    legend.text = element_text(size = 7, family = "mono")
  ) +
  labs(title = sprintf("(b) Family composition — assigned only (top %d + Other)", top_n))

tab_b2 <- raw %>%
  summarise(
    n_ass = sum(assigned_f),
    n_un = sum(!assigned_f)
  ) %>%
  pivot_longer(everything(), names_to = "grp", values_to = "n") %>%
  mutate(
    lab = ifelse(grp == "n_ass", "Assigned_F", "Unclassified_F"),
    frac = n / sum(n),
    pct = n / sum(n) * 100
  )

lv_b2 <- tab_b2$lab
lab_b2 <- legend_pct_count_name(as.character(tab_b2$lab), tab_b2$pct, tab_b2$n)

p_b2 <- ggplot(tab_b2, aes(x = "", y = frac, fill = lab)) +
  geom_col(width = 1, colour = "white", linewidth = 0.3) +
  coord_polar(theta = "y") +
  scale_fill_manual(
    values = c("Assigned_F" = "#6699CC", "Unclassified_F" = COL_UNASSIGNED),
    breaks = lv_b2,
    labels = lab_b2,
    name = NULL
  ) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 9.5, hjust = 0.5),
    legend.position = "bottom",
    legend.text = element_text(size = 7, family = "mono")
  ) +
  labs(title = "Assigned vs unclassified (family)")

p_b <- p_b1 / p_b2 + plot_layout(heights = c(1.1, 1))

# ---- (c) Chord: Top N class × Top N family (same collapse as (a)(b)) ----
mat_long <- raw %>%
  count(class_plot, family_plot, name = "value")

mat_long <- mat_long %>%
  mutate(
    from = paste0("C|", class_plot),
    to = paste0("F|", family_plot)
  )

sectors <- unique(c(mat_long$from, mat_long$to))
grid_col <- rep("#CCCCCC", length(sectors))
names(grid_col) <- sectors

cls_u <- unique(as.character(mat_long$class_plot))
fam_u <- unique(as.character(mat_long$family_plot))
lv_c <- levels(tab_a$class_plot)
lv_f <- levels(tab_f$family_plot)
pal_c_map <- assign_fig3_fill_colors(lv_c, "class")
pal_f_map <- assign_fig3_fill_colors(lv_f, "family")
pal_c <- pal_c_map[cls_u]
pal_c[is.na(pal_c)] <- COL_UNASSIGNED
names(pal_c) <- cls_u
pal_f_ch <- pal_f_map[fam_u]
pal_f_ch[is.na(pal_f_ch)] <- COL_UNASSIGNED
names(pal_f_ch) <- fam_u

for (cl in cls_u) {
  hits <- sectors[startsWith(sectors, paste0("C|", cl))]
  grid_col[hits] <- pal_c[[cl]]
}
for (fm in fam_u) {
  hits <- sectors[startsWith(sectors, paste0("F|", fm))]
  grid_col[hits] <- pal_f_ch[[fm]]
}

png_chord <- paste0(out_base, "_panel_c_chord.png")
pdf_chord <- paste0(out_base, "_panel_c_chord.pdf")

plot_chord_file <- function(path, w = 7.5, h = 7.5) {
  if (grepl("\\.pdf$", path)) {
    if (isTRUE(capabilities("cairo"))) {
      cairo_pdf(path, width = w, height = h)
    } else {
      pdf(path, width = w, height = h)
    }
  } else {
    png(path, width = w * 150, height = h * 150, res = 150)
  }
  par(mar = c(1, 1, 1, 1) * 0.1)
  circos.clear()
  chordDiagram(
    x = as.data.frame(mat_long %>% select(from, to, value)),
    grid.col = grid_col,
    transparency = 0.28,
    annotationTrack = "grid"
  )
  title(
    main = sprintf("(c) Class–family links — top %d + Other (vOTU counts)", top_n),
    cex.main = 1.05,
    font.main = 2
  )
  dev.off()
  circos.clear()
}

plot_chord_file(png_chord, 7.5, 7.5)
plot_chord_file(pdf_chord, 7.5, 7.5)

# Combine a + b into one PDF; c is a separate file (circlize is not ggplot)
combo <- p_a + p_b + plot_layout(widths = c(1.05, 1))
out_ab <- paste0(out_base, "_panels_a_b.pdf")
pdf_ab <- if (isTRUE(capabilities("cairo"))) cairo_pdf else pdf
ggsave(out_ab, combo, width = 11, height = 6.2, device = pdf_ab, dpi = 300)
ggsave(sub("\\.pdf$", ".png", out_ab), combo, width = 11, height = 6.2, dpi = 200, bg = "white")

cap <- c(
  "Fig.3-style panels from user TSV.",
  paste("Source:", normalizePath(tsv, winslash = "/", mustWork = FALSE)),
  paste("N rows used:", nrow(raw)),
  comp_note,
  "",
  sprintf("(a) Donut: class from taxonomy; show top %d + Other (Unclassified_C if class empty).", top_n),
  sprintf("(b) Top pie: assigned family only, top %d + Other; bottom: Assigned_F vs Unclassified_F.", top_n),
  "Legend (a)(b): monospace labels as  pct%%  count  name.",
  sprintf("(c) chordDiagram: same top %d + Other on class and family; sectors C| / F|.", top_n),
  "Colours: ranked palette + Unclassified grey + Other maroon (NM-style legend).",
  "Outputs:",
  basename(out_ab),
  basename(png_chord),
  basename(pdf_chord),
  "",
  "Note: NM-style outer tick scales on chord require extra circos.track; extend as needed."
)
writeLines(cap, paste0(out_base, "_caption.txt"))

message("[OK] ", normalizePath(out_ab, winslash = "/"))
message("[OK] ", normalizePath(png_chord, winslash = "/"))
