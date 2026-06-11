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
# Top-N representative vOTUs per enterotype; block heatmap
#
# Marker scoring (--marker-mode):
#   vs_max_other (default): cluster-specific contrast × prevalence bonus; absolute abundance not emphasized.
#   composite: vs_max_other × pmax(log10(mu_k+p)+6, 0.05) for brighter within-cluster signal.
#   fold_max_other: log10((mu_k+p)/(max_other+p)) × prevalence bonus.
#   delta_mean_rest: mu_k - mean(non-k samples) (legacy).
# Reasonable filters (combine to avoid flat or all-black rows):
#   1) Candidate pool: tune --prev / --top-taxa like main enterotype run.
#   2) Cluster-specific: composite or vs_max_other + --max-mean-any-other.
#   3) Absolute brightness: --min-mean-in-et; --min-prev-in-et.
#   4) Plot columns: --samples-per-et / --max-columns to avoid overcrowding.
# Heatmap: rows blocked by enterotype; cols = samples sorted enterotype → sample_id; color log10(rel+pseudo).
# Also: {out_prefix}_marker_abundance_boxplot.pdf/png — facets by enterotype.
#
# Reduce figure size:
#   - Rows: --top-n (default 10) markers per enterotype.
#   - Columns = samples in heatmap. Either:
#       --samples-per-et 40   up to 40 random samples per ET;
#       --max-columns 200     cap total columns at 200 across ETs.
#
# Usage (Phage root):
#   Rscript meta_qc/plot_enterotype_marker_heatmap.R \
#     --assign meta_qc/Enterotype_vOTU/PAM_k2_samples.tsv
#
#   Rscript meta_qc/plot_enterotype_marker_heatmap.R \
#     --assign meta_qc/Enterotype_vOTU/PAM_k4_samples.tsv \
#     --marker-mode vs_max_other --top-n 10 --samples-per-et 35 \
#     --min-depth 3000 --prev 0.05 --top-taxa 8000 \
#     --out-prefix meta_qc/Enterotype_vOTU/marker_hm_k4
#
#   Brightness + specificity: --marker-mode composite --min-mean-in-et 1e-5 --min-prev-in-et 0.15
#   Stricter other-ET dim: --max-mean-any-other 3e-5
#   Suppress pan-core bright vOTUs: --min-fold-k-vs-max-other 1.5
#     Require mean_rel_k >= fold × (max(other ET means)+pseudo).
#
# Overlay taxonomy / host / virus_class on heatmap row labels:
#   1) Prepare TSV aligned to votu (e.g. meta_qc/marker_list_20_by_enterotype_taxonomy_host.tsv
#      or annotate_marker_votu_taxonomy_host.R output); columns must include votu;
#      optional: best_host_genus, virus_class, virus_family_rank, taxonomy, best_host_taxonomy.
#   2) Run with:
#        --annot-tsv path/to/annot.tsv --heatmap-ylabels annot
#      Y-axis two lines: ETk + truncated vOTU id; (host genus | virus_class).
#      Unmatched rows use "-" / "?".
#
# Requires: data.table, ggplot2
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
    out_prefix = getv("--out-prefix", ""),
    top_n = as.integer(getv("--top-n", "10")),
    min_depth = as.integer(getv("--min-depth", "3000")),
    prev = as.numeric(getv("--prev", "0.05")),
    top_taxa = as.integer(getv("--top-taxa", "5000")),
    pseudo = as.numeric(getv("--pseudo", "1e-7")),
    samples_per_et = {
      s <- getv("--samples-per-et", "")
      if (!nzchar(s)) NA_integer_ else as.integer(s)
    },
    max_samples_plot = {
      s1 <- getv("--max-samples-plot", "")
      s2 <- getv("--max-columns", "")
      s <- if (nzchar(s2)) s2 else s1
      if (!nzchar(s)) NA_integer_ else as.integer(s)
    },
    seed = as.integer(getv("--seed", "123")),
    marker_mode = tolower(trimws(getv("--marker-mode", "vs_max_other"))),
    min_mean_in_et = as.numeric(getv("--min-mean-in-et", "0")),
    max_mean_any_other = as.numeric(getv("--max-mean-any-other", "Inf")),
    min_prev_in_et = {
      s <- getv("--min-prev-in-et", "")
      if (!nzchar(s)) NA_real_ else as.numeric(s)
    },
    min_fold_k_vs_max_other = {
      s <- getv("--min-fold-k-vs-max-other", "")
      if (!nzchar(s)) NA_real_ else as.numeric(s)
    },
    annot_tsv = getv("--annot-tsv", ""),
    heatmap_ylabels = tolower(trimws(getv("--heatmap-ylabels", "id_only")))
  )
}

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
root <- if (length(script_path)) dirname(dirname(normalizePath(script_path))) else getwd()
if (!file.exists(file.path(root, "pvca_Y_count.tsv"))) root <- getwd()

pa <- parse_args()

assign_path <- if (file.exists(pa$assign_path)) normalizePath(pa$assign_path) else file.path(root, pa$assign_path)

count_path <- if (nzchar(pa$count_path)) {
  if (file.exists(pa$count_path)) normalizePath(pa$count_path) else file.path(root, pa$count_path)
} else {
  file.path(root, "pvca_Y_count.tsv")
}

out_prefix <- if (nzchar(pa$out_prefix)) {
  normalizePath(pa$out_prefix, winslash = "/", mustWork = FALSE)
} else {
  file.path(dirname(assign_path), "enterotype_marker_heatmap")
}
dir.create(dirname(out_prefix), showWarnings = FALSE, recursive = TRUE)

# Enterotype top bar colors: match PCoA (ET1 dark blue, ET2 light blue, ET3 dark orange, ET4 light orange)
ET_ANN_PAL <- c(
  "#1f78b4", "#a6cee3", "#ff7f0e", "#fdbf6f",
  "#33a02c", "#b15928", "#6a3d9a", "#cab2d6", "#fb9a99", "#a6cee3"
)
# Heatmap cells: log10 rel. abundance low→blue, mid→white, high→orange
HM_FILL_COLS <- c("#2166ac", "#6baed6", "#c6dbef", "#f7f7f7", "#fee8c8", "#fdae6b", "#fd8d3c", "#e6550d")

asg[, sample_id := as.character(sample_id)]
asg[, cluster := as.integer(cluster)]

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
rm(cnt)

depth <- rowSums(mat)
keep <- depth >= pa$min_depth
mat <- mat[keep, , drop = FALSE]
asg <- asg[match(rownames(mat), asg$sample_id)]

prev <- colSums(mat > 0) / nrow(mat)
mat <- mat[, prev >= pa$prev, drop = FALSE]

rel0 <- mat / rowSums(mat)
mean_rel <- colMeans(rel0)
if (is.finite(pa$top_taxa) && !is.na(pa$top_taxa) && pa$top_taxa > 0L && ncol(mat) > pa$top_taxa) {
  ord <- order(mean_rel, decreasing = TRUE)
  take <- ord[seq_len(pa$top_taxa)]
  mat <- mat[, take, drop = FALSE]
  rel0 <- rel0[, take, drop = FALSE]
}

ks <- sort(unique(asg$cluster))

mu_mat <- sapply(ks, function(ll) colMeans(rel0[asg$cluster == ll, , drop = FALSE]))
rownames(mu_mat) <- colnames(rel0)
colnames(mu_mat) <- as.character(ks)
prev_mat <- sapply(ks, function(ll) colMeans(rel0[asg$cluster == ll, , drop = FALSE] > 0))
rownames(prev_mat) <- colnames(rel0)
colnames(prev_mat) <- as.character(ks)

mode <- pa$marker_mode

min_et <- pa$min_mean_in_et
if (is.na(min_et)) min_et <- 0
max_other_th <- pa$max_mean_any_other
if (is.na(max_other_th)) max_other_th <- Inf
min_prev_et <- pa$min_prev_in_et
if (length(min_prev_et) != 1L || is.na(min_prev_et)) min_prev_et <- NA_real_
min_fold_k <- pa$min_fold_k_vs_max_other
if (length(min_fold_k) != 1L || is.na(min_fold_k)) min_fold_k <- NA_real_

topn <- pa$top_n
if (is.na(topn) || topn < 1L) topn <- 10L

marker_rows <- list()
score_tbl <- list()

for (k in ks) {
  ik <- asg$cluster == k
  n_k <- sum(ik)
  n_o <- sum(!ik)
  if (n_k < 3L || n_o < 3L) {
    message("Waring: cluster ", k, " sample too litter (k=", n_k, ", rest=", n_o, ")")
  }
  sk <- as.character(k)
  others <- setdiff(as.character(ks), sk)
  if (length(others) >= 1L) {
    mu_max_other <- apply(mu_mat[, others, drop = FALSE], 1L, max)
    prev_max_other <- apply(prev_mat[, others, drop = FALSE], 1L, max)
  } else {
    mu_max_other <- rep(0, nrow(mu_mat))
    prev_max_other <- rep(0, nrow(prev_mat))
  }
  names(mu_max_other) <- rownames(mu_mat)
  names(prev_max_other) <- rownames(prev_mat)

  Rk <- rel0[ik, , drop = FALSE]
  Ro <- rel0[!ik, , drop = FALSE]
  mu_o_out <- colMeans(Ro)
  po_out <- colMeans(Ro > 0)

  if (mode == "delta_mean_rest") {
    mu_k <- colMeans(Rk)
    pk <- colMeans(Rk > 0)
    po <- po_out
    score <- (mu_k - mu_o_out) * (1 + pmax(pk - po, 0))
    names(score) <- colnames(rel0)
  } else {
    mu_k <- mu_mat[, sk]
    pk <- prev_mat[, sk]
    names(mu_k) <- rownames(mu_mat)
    names(pk) <- rownames(prev_mat)
    pm <- pa$pseudo
    prev_bonus <- 1 + pmax(pk - prev_max_other, 0)
    base_contrast <- (mu_k - mu_max_other) * prev_bonus
    if (mode == "vs_max_other") {
      score <- base_contrast
    } else if (mode == "composite") {
      # Boost within-cluster absolute abundance weight to avoid all-black heatmap
      score <- base_contrast * pmax(log10(mu_k + pm) + 6, 0.05)
    } else {
      score <- log10((mu_k + pm) / (mu_max_other + pm)) * prev_bonus
    }
    names(score) <- names(mu_k)
  }

  ok <- is.finite(score) & is.finite(mu_k) & (mu_k >= min_et)
  if (is.finite(min_prev_et) && min_prev_et > 0) {
    ok <- ok & is.finite(pk) & (pk >= min_prev_et)
  }
  if (is.finite(min_fold_k) && min_fold_k > 0 && length(others) >= 1L) {
    ok <- ok & is.finite(mu_max_other) & (mu_k >= min_fold_k * (mu_max_other + pa$pseudo))
  }
  if (is.finite(max_other_th)) {
    ok <- ok & is.finite(mu_max_other) & (mu_max_other <= max_other_th)
  }
  o_all <- order(score, decreasing = TRUE)
  o_f <- o_all[ok[o_all]]
  if (length(o_f) < topn) {
    message()
    need <- topn - length(o_f)
    rest <- setdiff(o_all, o_f)
    o_take <- c(o_f, rest[seq_len(min(need, length(rest)))])
  } else {
    o_take <- o_f[seq_len(topn)]
  }
  take <- o_take[seq_len(min(topn, length(o_take)))]
  votu_k <- names(score)[take]
  mu_k_v <- as.numeric(mu_k[votu_k])
  pk_v <- as.numeric(pk[votu_k])
  mmx_v <- as.numeric(mu_max_other[votu_k])
  pmx_v <- as.numeric(prev_max_other[votu_k])
  muo_v <- as.numeric(mu_o_out[votu_k])
  po_v <- as.numeric(po_out[votu_k])
  marker_rows[[length(marker_rows) + 1L]] <- data.table(
    enterotype = k,
    rank = seq_along(votu_k),
    votu = votu_k,
    score = as.numeric(score[votu_k]),
    mean_rel_k = mu_k_v,
    mean_rel_max_other_et = mmx_v,
    mean_rel_pooled_rest = muo_v,
    prev_k = pk_v,
    prev_max_other_et = pmx_v,
    prev_pooled_rest = po_v,
    marker_mode = mode
  )
  score_tbl[[length(score_tbl) + 1L]] <- data.table(
    enterotype = k,
    votu = names(score),
    score = as.numeric(score),
    mean_rel_k = as.numeric(mu_k[names(score)]),
    mean_rel_max_other_et = as.numeric(mu_max_other[names(score)]),
    mean_rel_pooled_rest = as.numeric(mu_o_out[names(score)]),
    prev_k = as.numeric(pk[names(score)]),
    prev_max_other_et = as.numeric(prev_max_other[names(score)]),
    prev_pooled_rest = as.numeric(po_out[names(score)]),
    marker_mode = mode
  )
}

markers_dt <- rbindlist(marker_rows)
if (nzchar(pa$annot_tsv)) {
  ap <- if (file.exists(pa$annot_tsv)) normalizePath(pa$annot_tsv) else file.path(root, pa$annot_tsv)
  if (file.exists(ap)) {
    ann_m <- fread(ap, sep = "\t", encoding = "UTF-8")
    if ("votu" %in% names(ann_m)) {
      ann_m[, votu := as.character(votu)]
      markers_dt[, votu := as.character(votu)]
      markers_dt[, `.marker_row_ord` := seq_len(.N)]
      add_cols <- intersect(
        c("taxonomy", "virus_class", "virus_family_rank", "best_host_genus", "best_host_taxonomy"),
        names(ann_m)
      )
      ann_u <- unique(ann_m[, c("votu", add_cols), with = FALSE], by = "votu")
      markers_dt <- merge(markers_dt, ann_u, by = "votu", all.x = TRUE, sort = FALSE)
      setorder(markers_dt, `.marker_row_ord`)
      markers_dt[, `.marker_row_ord` := NULL]
    } else {
      message()
    }
  } else {
    message()
  }
}
fwrite(markers_dt, paste0(out_prefix, "_top_markers.tsv"), sep = "\t")
fwrite(rbindlist(score_tbl), paste0(out_prefix, "_all_scores_long.tsv"), sep = "\t")


# ---- Marker vOTU rel. abundance boxplot: facet by enterotype ----
bx_rows <- list()
for (i in seq_len(nrow(markers_dt))) {
  et <- markers_dt$enterotype[i]
  v <- markers_dt$votu[i]
  rk <- markers_dt$rank[i]
  if (!v %in% colnames(rel0)) next
  sid_et <- asg[cluster == et, sample_id]
  sid_et <- sid_et[sid_et %in% rownames(rel0)]
  if (length(sid_et) < 2L) next
  val <- as.numeric(rel0[sid_et, v, drop = TRUE])
  bx_rows[[length(bx_rows) + 1L]] <- data.table(
    enterotype = et,
    votu = v,
    rank = rk,
    rel = val
  )
}
bx <- rbindlist(bx_rows)
if (nrow(bx) > 0L) {
  bx[, enterotype := factor(enterotype, levels = ks)]
  for (k2 in ks) {
    lv <- markers_dt[enterotype == k2][order(rank)]$votu
    if (length(lv)) bx[enterotype == k2, votu_f := factor(votu, levels = lv)]
  }
  p_box <- ggplot(bx, aes(x = votu_f, y = rel)) +
    geom_boxplot(
      width = 0.68,
      outlier.size = 0.45,
      linewidth = 0.35,
      fill = "#eef5fc",
      colour = "#2166ac",
      alpha = 0.95
    ) +
    facet_wrap(~enterotype, scales = "free_x", ncol = min(2L, length(ks))) +
    labs(
      title = "Marker vOTU relative abundance (within enterotype samples)",
      x = "vOTU (markers ordered by rank within each ET)",
      y = "Relative abundance"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 50, hjust = 1, size = 6.2, vjust = 1),
      strip.text = element_text(face = "bold", size = 11),
      strip.background = element_rect(fill = "#e8f2fa", colour = NA),
      plot.title = element_text(size = 12, hjust = 0.5),
      panel.grid.major.x = element_blank()
    )
  n_et <- length(ks)
  w_box <- min(20, max(9.5, pa$top_n * 0.72 * max(2L, ceiling(n_et / 2))))
  h_box <- max(5, 3.4 + ceiling(n_et / 2) * 2.35)
  ggsave(paste0(out_prefix, "_marker_abundance_boxplot.pdf"), p_box, width = w_box, height = h_box, limitsize = FALSE)
  ggsave(paste0(out_prefix, "_marker_abundance_boxplot.png"), p_box, width = w_box, height = h_box, dpi = 300, limitsize = FALSE)
} else {
  message()
}

markers_dt[, block := factor(paste0("ET", enterotype), levels = paste0("ET", ks))]
yl_mode <- pa$heatmap_ylabels
if (!yl_mode %in% c("id_only", "annot")) {
  message("Unknown --heatmap-ylabels=", yl_mode, "，change to id_only")
  yl_mode <- "id_only"
}
if (yl_mode == "annot") {
  if (!"best_host_genus" %in% names(markers_dt)) markers_dt[, best_host_genus := NA_character_]
  if (!"virus_class" %in% names(markers_dt)) markers_dt[, virus_class := NA_character_]
  trunc_str <- function(s, n = 40L) {
    s <- as.character(s)
    ok <- !is.na(s) & nzchar(trimws(s))
    out <- s
    w <- ok & nchar(trimws(s)) > n
    t <- trimws(s)
    out[w] <- paste0(substr(t[w], 1L, max(1L, n - 3L)), "...")
    out[!ok] <- NA_character_
    out
  }
  markers_dt[, ylabel := {
    etl <- paste0("ET", enterotype)
    sv <- trunc_str(votu, 42L)
    hg <- trimws(as.character(best_host_genus))
    hg[!nzchar(hg) | is.na(hg)] <- "-"
    vc <- trimws(as.character(virus_class))
    vc[!nzchar(vc) | is.na(vc)] <- "?"
    vc <- trunc_str(vc, 28L)
    paste0(etl, ": ", sv, "\n(", hg, " | ", vc, ")")
  }]
} else {
  markers_dt[, ylabel := paste0("ET", enterotype, ": ", votu)]
}
lev_y <- markers_dt$ylabel
markers_dt[, y := factor(ylabel, levels = lev_y)]
use_annot_y <- yl_mode == "annot"

asg[, et_col := factor(cluster, levels = ks)]
sid_ord <- asg$sample_id[order(asg$cluster, asg$sample_id)]
set.seed(pa$seed)
plot_sub <- is.finite(pa$samples_per_et) && !is.na(pa$samples_per_et) && pa$samples_per_et > 0L
tot_sub <- is.finite(pa$max_samples_plot) && !is.na(pa$max_samples_plot) && pa$max_samples_plot > 0L

if (isTRUE(plot_sub)) {
  cap <- pa$samples_per_et
  sid_keep <- unlist(lapply(ks, function(kk) {
    s <- asg[cluster == kk, sample_id]
    if (length(s) <= cap) return(s)
    sample(s, cap)
  }))
  sid_ord <- sid_ord[sid_ord %in% sid_keep]
  rel_plot <- rel0[sid_ord, , drop = FALSE]
  asg_plot <- asg[match(sid_ord, sample_id)]
} else if (isTRUE(tot_sub) && length(sid_ord) > pa$max_samples_plot) {
  per <- max(1L, floor(pa$max_samples_plot / length(ks)))
  sid_keep <- unlist(lapply(ks, function(kk) {
    s <- asg[cluster == kk, sample_id]
    if (length(s) <= per) return(s)
    sample(s, per)
  }))
  sid_ord <- sid_ord[sid_ord %in% sid_keep]
  rel_plot <- rel0[sid_ord, , drop = FALSE]
  asg_plot <- asg[match(sid_ord, sample_id)]
} else {
  rel_plot <- rel0[sid_ord, , drop = FALSE]
  asg_plot <- asg[match(sid_ord, sample_id)]
}

# After subsample, keep enterotype → sample_id column order
asg_plot <- asg_plot[order(asg_plot$cluster, asg_plot$sample_id)]
sid_ord <- asg_plot$sample_id
rel_plot <- rel0[sid_ord, , drop = FALSE]

# Expand per marker row to avoid merge inflation when vOTU in multiple ET tops
ml <- rbindlist(lapply(seq_len(nrow(markers_dt)), function(i) {
  v <- markers_dt$votu[i]
  data.table(
    sample_id = sid_ord,
    votu = v,
    ylabel = markers_dt$ylabel[i],
    block = markers_dt$block[i],
    enterotype = markers_dt$enterotype[i],
    log10_rel = log10(rel_plot[, v] + pa$pseudo)
  )
}))
ml[, sample_id := factor(sample_id, levels = sid_ord)]
ml[, y := factor(ylabel, levels = lev_y)]

ann <- unique(asg_plot[, .(sample_id, cluster)])
ann[, sample_id := factor(sample_id, levels = levels(ml$sample_id))]
uk_ann <- sort(unique(as.integer(ann$cluster)))
pal_ann <- ET_ANN_PAL[((uk_ann - 1L) %% length(ET_ANN_PAL)) + 1L]
names(pal_ann) <- as.character(uk_ann)

# Column annotation bar 1: enterotype
p_ann_et <- ggplot(ann, aes(x = sample_id, y = 1L, fill = factor(cluster))) +
  geom_tile() +
  scale_fill_manual(values = pal_ann, breaks = names(pal_ann), name = "Enterotype") +
  scale_y_continuous(expand = c(0, 0), breaks = NULL, labels = NULL) +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(0, 5.5, 0, 5.5))

br_fill <- as.numeric(quantile(ml$log10_rel, probs = c(0.02, 0.98), na.rm = TRUE))
if (!is.finite(br_fill[1])) br_fill[1] <- min(ml$log10_rel, na.rm = TRUE)
if (!is.finite(br_fill[2])) br_fill[2] <- max(ml$log10_rel, na.rm = TRUE)
if (br_fill[2] <= br_fill[1]) br_fill[2] <- br_fill[1] + 1e-6

p_hm <- ggplot(ml, aes(x = sample_id, y = y, fill = log10_rel)) +
  geom_raster() +
  facet_grid(block ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradientn(
    colours = HM_FILL_COLS,
    limits = br_fill,
    oob = scales::squish,
    name = "log10(rel+p)",
    na.value = "#eeeeee"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(
    x = NULL,
    y = NULL,
    title = "Enterotype markers (columns: ET, then sample)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(
      size = if (isTRUE(use_annot_y)) 4.6 else 5.5,
      lineheight = if (isTRUE(use_annot_y)) 0.95 else 1
    ),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "#e8f2fa", colour = NA),
    panel.spacing.y = unit(0.35, "lines"),
    plot.title = element_text(hjust = 0.5, size = 12)
  )

# Combine: enterotype top bar + heatmap
if (!requireNamespace("patchwork", quietly = TRUE)) {
  n_s <- length(levels(ml$sample_id))
  n_r <- length(levels(ml$y))
  row_h <- if (isTRUE(use_annot_y)) 0.32 else 0.22
  h_main <- max(5, n_r * row_h)
  w_main <- min(if (isTRUE(use_annot_y)) 20 else 18, max(6.5, n_s * 0.022))
  ggsave(paste0(out_prefix, "_heatmap.pdf"), p_hm, width = w_main, height = h_main, limitsize = FALSE)
  ggsave(paste0(out_prefix, "_heatmap.png"), p_hm, width = w_main, height = h_main, dpi = 220, limitsize = FALSE)
} else {
  n_s <- length(levels(ml$sample_id))
  n_r <- length(levels(ml$y))
  h_ann <- 0.45
  row_h <- if (isTRUE(use_annot_y)) 0.32 else 0.22
  h_main <- max(5, n_r * row_h)
  w_main <- min(if (isTRUE(use_annot_y)) 20 else 18, max(6.5, n_s * 0.022))
  pw <- patchwork::wrap_plots(p_ann_et, p_hm, ncol = 1, heights = c(h_ann, h_main))
  h_tot <- h_main + h_ann + 0.3
  ggsave(paste0(out_prefix, "_heatmap.pdf"), pw, width = w_main, height = h_tot, limitsize = FALSE)
  ggsave(paste0(out_prefix, "_heatmap.png"), pw, width = w_main, height = h_tot, dpi = 220, limitsize = FALSE)
}