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
# By enterotype: ratio of vOTUs with Lactobacillus/Limosilactobacillus vs Prevotella hosts
# (Prevotella, exact best_host_genus match) relative abundance sums
#
# Per sample:
#   lact_rel  = sum_j rel_ij, j with best_host_genus in {Lactobacillus,Limosilactobacillus}
#   prev_rel  = sum_j rel_ij, j with best_host_genus == Prevotella
#   log2_ratio = log2((lact_rel + pseudo) / (prev_rel + pseudo))
#
# Plot: X=enterotype, Y=log2_ratio; violin + narrow box (alpha diversity style, no jitter).
#   Annotate KW global p; pairwise Wilcoxon with BH (--no-pairwise-p for global only).
#
# Input:
#   - PAM_k*_samples.tsv: sample_id, cluster
#   - pvca_Y_count.tsv: count matrix
#   - Annotation table (default meta_qc/caudoviricetes_votu_annotated_merged.tsv): vOTU, best_host_genus
#
# Usage (Phage root):
#   Rscript meta_qc/plot_enterotype_lab_prev_ratio_boxplot.R \
#     --assign meta_qc/Enterotype_vOTU_hellinger/PAM_k4_samples.tsv \
#     --out-prefix meta_qc/Enterotype_vOTU_hellinger/lab_prev_ratio_k4
#
#   Rscript meta_qc/plot_enterotype_lab_prev_ratio_boxplot.R \
#     --assign meta_qc/Enterotype_vOTU/PAM_k4_samples.tsv \
#     --annot-tsv meta_qc/caudoviricetes_votu_annotated_merged.tsv \
#     --min-depth 3000 --pseudo 1e-8
#
#   Global Kruskal–Wallis only, no pairwise brackets:
#     ... --no-pairwise-p
#
# Requires: data.table, ggplot2, ggpubr (KW + pairwise Wilcoxon, BH)
#
# Reproducibility check (if plots differ, compare these files):
#   {out_prefix}_plot_run_manifest.tsv   paths, MD5, min_depth, pseudo, seed
#   {out_prefix}_plot_summary_by_enterotype.tsv  n and log2 median/mean per ET
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

parse_args <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  getv <- function(flag, default = NA_character_) {
    i <- match(flag, a)
    if (is.na(i) || i >= length(a)) return(default)
    a[[i + 1]]
  }
  list(
    assign_path = getv("--assign", ""),
    count_path = getv("--count-path", ""),
    annot_path = getv("--annot-tsv", ""),
    out_prefix = getv("--out-prefix", ""),
    min_depth = as.integer(getv("--min-depth", "0")),
    pseudo = as.numeric(getv("--pseudo", "1e-8")),
    plot_format = tolower(getv("--plot-format", "pdf")),
    seed = as.integer(getv("--seed", "123")),
    pairwise_p = {
      a <- commandArgs(trailingOnly = TRUE)
      !("--no-pairwise-p" %in% a)
    }
  )
}

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
root <- if (length(script_path)) dirname(dirname(normalizePath(script_path))) else getwd()
root_inferred <- root

pa <- parse_args()

assign_path <- if (file.exists(pa$assign_path)) normalizePath(pa$assign_path) else file.path(root, pa$assign_path)


count_path <- if (nzchar(pa$count_path)) {
  if (file.exists(pa$count_path)) normalizePath(pa$count_path) else file.path(root, pa$count_path)
} else {
  file.path(root, "pvca_Y_count.tsv")
}

annot_path <- if (nzchar(pa$annot_path)) {
  if (file.exists(pa$annot_path)) normalizePath(pa$annot_path) else file.path(root, pa$annot_path)
} else {
  file.path(root, "meta_qc", "caudoviricetes_votu_annotated_merged.tsv")
}


out_prefix <- if (nzchar(pa$out_prefix)) {
  normalizePath(pa$out_prefix, winslash = "/", mustWork = FALSE)
} else {
  file.path(dirname(assign_path), "lab_prev_ratio_by_enterotype")
}
dir.create(dirname(out_prefix), showWarnings = FALSE, recursive = TRUE)

ET_PAL <- c(
  "#1f78b4", "#a6cee3", "#ff7f0e", "#fdbf6f",
  "#33a02c", "#b15928", "#6a3d9a", "#cab2d6", "#fb9a99", "#a6cee3"
)

LAB_GENERA <- c("Lactobacillus", "Limosilactobacillus")
PREV_GENUS <- "Prevotella"


an <- fread(annot_path, sep = "\t", encoding = "UTF-8")
vcol <- if ("vOTU" %in% names(an)) "vOTU" else if ("votu" %in% names(an)) "votu" else stop()

an[, votu := as.character(get(vcol))]
an[, hg := trimws(as.character(best_host_genus))]
an <- an[nzchar(votu)]

lab_v <- unique(an[hg %in% LAB_GENERA, votu])
prev_v <- unique(an[hg == PREV_GENUS, votu])

asg <- fread(assign_path, sep = "\t", encoding = "UTF-8")
asg[, sample_id := as.character(sample_id)]
asg[, cluster := as.integer(cluster)]
ks <- sort(unique(asg$cluster))
asg[, et_label := factor(paste0("ET", cluster), levels = paste0("ET", ks))]


cnt <- fread(count_path, sep = "\t", encoding = "UTF-8", showProgress = TRUE)

cnt[, sample_id := as.character(sample_id)]
cnt <- cnt[sample_id %in% asg$sample_id]
asg <- asg[sample_id %in% cnt$sample_id]
setkey(asg, sample_id)
asg <- asg[J(cnt$sample_id)]
stopifnot(all(asg$sample_id == cnt$sample_id))

feat_cols <- setdiff(names(cnt), "sample_id")
mat <- as.matrix(cnt[, ..feat_cols])
storage.mode(mat) <- "double"
rownames(mat) <- cnt$sample_id
colnames(mat) <- feat_cols

depth <- rowSums(mat)
if (is.finite(pa$min_depth) && !is.na(pa$min_depth) && pa$min_depth > 0L) {
  keep <- depth >= pa$min_depth
  mat <- mat[keep, , drop = FALSE]
  asg <- asg[match(rownames(mat), sample_id), ]
}

rel <- mat / rowSums(mat)
lab_in <- intersect(lab_v, colnames(rel))
prev_in <- intersect(prev_v, colnames(rel))

lact_rel <- rowSums(rel[, lab_in, drop = FALSE])
prev_rel <- rowSums(rel[, prev_in, drop = FALSE])
eps <- pa$pseudo
if (!is.finite(eps) || eps <= 0) eps <- 1e-8
log2_ratio <- log2((lact_rel + eps) / (prev_rel + eps))

out_dt <- data.table(
  sample_id = rownames(rel),
  enterotype = as.integer(asg$cluster[match(rownames(rel), asg$sample_id)]),
  et_label = asg$et_label[match(rownames(rel), asg$sample_id)],
  depth = rowSums(mat),
  lact_host_rel = lact_rel,
  prev_host_rel = prev_rel,
  n_lab_votu_detected = rowSums(mat[, lab_in, drop = FALSE] > 0),
  n_prev_votu_detected = rowSums(mat[, prev_in, drop = FALSE] > 0),
  log2_lact_prev_ratio = log2_ratio
)
fwrite(out_dt, paste0(out_prefix, "_per_sample_lab_prev.tsv"), sep = "\t")

# Run manifest: check whether plot differences come from input paths/file versions (p-value annotations do not change these stats)
sum_by_et <- out_dt[, .(
  n = .N,
  median_log2 = median(log2_lact_prev_ratio, na.rm = TRUE),
  mean_log2 = mean(log2_lact_prev_ratio, na.rm = TRUE)
), by = et_label][order(et_label)]
prov <- data.table(
  field = c(
    "timestamp", "root_used", "assign_path", "count_path", "annot_path",
    "min_depth", "pseudo", "seed", "pairwise_p",
    "md5_assign", "md5_count", "md5_annot"
  ),
  value = c(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    root,
    assign_path,
    count_path,
    annot_path,
    as.character(pa$min_depth),
    as.character(pa$pseudo),
    as.character(pa$seed),
    as.character(isTRUE(pa$pairwise_p)),
    as.character(tools::md5sum(assign_path)),
    as.character(tools::md5sum(count_path)),
    as.character(tools::md5sum(annot_path))
  )
)
fwrite(prov, paste0(out_prefix, "_plot_run_manifest.tsv"), sep = "\t", quote = FALSE)
fwrite(sum_by_et, paste0(out_prefix, "_plot_summary_by_enterotype.tsv"), sep = "\t", quote = FALSE)

pal_use <- ET_PAL[((ks - 1L) %% length(ET_PAL)) + 1L]
names(pal_use) <- paste0("ET", ks)

tab_n <- as.integer(table(factor(out_dt$enterotype, levels = ks)))
cap_sz <- paste(paste0("ET", ks, " n=", tab_n), collapse = "; ")

set.seed(pa$seed)
rng <- range(out_dt$log2_lact_prev_ratio, na.rm = TRUE)
y_span <- diff(rng)
if (!is.finite(y_span) || y_span < 1e-12) y_span <- 1
et_lv <- levels(out_dt$et_label)
n_pairs <- if (length(et_lv) < 2L) 0L else as.integer(ncol(combn(et_lv, 2L)))
# Reserve y-axis space above for global and pairwise brackets
top_pad <- if (length(ks) < 2L) {
  0.1
} else if (!isTRUE(pa$pairwise_p) || n_pairs <= 1L) {
  0.14
} else if (n_pairs <= 3L) {
  0.22
} else {
  min(0.38, 0.16 + 0.028 * n_pairs)
}

p <- ggplot(out_dt, aes(x = et_label, y = log2_lact_prev_ratio, fill = et_label)) +
  geom_violin(trim = FALSE, alpha = 0.55, colour = NA, width = 0.85) +
  geom_boxplot(width = 0.12, outlier.size = 0.35, linewidth = 0.35, alpha = 0.92, colour = "grey25") +
  scale_fill_manual(values = pal_use, name = "Enterotype") +
  labs(
    title = "Lactobacillus/Limosilactobacillus vs Prevotella host-associated vOTUs",
    subtitle = paste0(
      "Per-sample sum of rel.abund; log2((Lact+Limosi+eps)/(Prevotella+eps)); eps=", format(eps, scientific = TRUE)
    ),
    caption = paste0(
      "Sample sizes: ", cap_sz,
      if (length(ks) >= 2L) {
        ". Stats: Kruskal–Wallis (global); pairwise Wilcoxon rank-sum, BH-adjusted p (two-sided)."
      } else {
        "."
      }
    ),
    x = NULL,
    y = "log2[(sum rel. Lactobacillus+Limosilactobacillus host vOTUs) / (sum rel. Prevotella host vOTUs)]"
  ) +
  scale_y_continuous(expand = ggplot2::expansion(mult = c(0.06, top_pad))) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 8.5, colour = "grey30", lineheight = 1.12),
    plot.caption = element_text(hjust = 0.5, size = 7.2, colour = "grey25"),
    plot.margin = ggplot2::margin(5.5, 8, 5.5, 8)
  )

if (length(ks) >= 2L) {
  if (isTRUE(pa$pairwise_p) && n_pairs >= 1L) {
    cm <- combn(et_lv, 2L)
    my_cmp <- lapply(seq_len(ncol(cm)), function(j) as.character(cm[, j]))
    p <- p +
      ggpubr::stat_compare_means(
        comparisons = my_cmp,
        method = "wilcox.test",
        p.adjust.method = "BH",
        label = "p.format",
        size = 3,
        bracket.size = 0.28,
        tip.length = 0.01
      )
  }
  kw_y <- rng[2] + y_span * (if (isTRUE(pa$pairwise_p) && n_pairs >= 1L) {
    0.11 + 0.036 * n_pairs
  } else {
    0.08
  })
  p <- p + ggpubr::stat_compare_means(method = "kruskal.test", label.y = kw_y)
}

fmt <- pa$plot_format
if (!fmt %in% c("pdf", "png")) fmt <- "pdf"
w <- min(8.5, max(5.2, 1.1 * length(ks)))
plot_h <- 5.2 + if (length(ks) <= 2L) 0 else if (length(ks) == 3L) 0.35 else 0.85
if (fmt == "pdf") {
  ggsave(paste0(out_prefix, "_boxplot_log2_ratio.pdf"), p, width = w, height = plot_h, limitsize = FALSE)
} else {
  ggsave(paste0(out_prefix, "_boxplot_log2_ratio.png"), p, width = w, height = plot_h, dpi = 300, limitsize = FALSE)
}