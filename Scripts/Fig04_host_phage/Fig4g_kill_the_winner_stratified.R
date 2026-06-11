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
# Species TPM vs virulent phage abundance infecting that species: per-sample match, per-species correlation + multiple testing
# Species abundance: taxonomy.species.tpm.txt (columns bioproject_runid, underscore-separated)
# Phage abundance: FINAL_all_projects_TPM_matrix.tsv (columns runid_TPM; strip _TPM for runid)
# Match columns by runid; Spearman per species; FDR correction
# Part 2: recompute correlation from per_sample for sig. negative species (signal samples) + forest plot + phage-stratified bars
# =============================================================================

library(data.table)
library(ggplot2)
if (requireNamespace("gridExtra", quietly = TRUE)) library(gridExtra)

# Working directory
if (dir.exists("D:\\F\\MicrobiomeMeta\\global\\Figure\\Phage")) # setwd: using PHAGE_ROOT above
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path)) tryCatch(setwd(dirname(script_path)), error = function(e) NULL)

# -----------------------------------------------------------------------------
# Input files
# -----------------------------------------------------------------------------
file_species_tpm   <- "taxonomy.species.tpm.txt"
file_phage_tpm     <- "FINAL_all_projects_TPM_matrix.tsv"
file_host_list     <- "host_infection_vOTU_list_per_host_species.csv"
file_bacphlip      <- "all_vOTUs_lifestyle_predictions.tsv"
file_phatyp        <- "high_quality_combineFilter.phatyp.csv"

strip_taxon_prefix <- function(x) sub("^[a-z]{1,2}__", "", trimws(as.character(x)))

# runids to exclude (not used in correlation or plots)
exclude_runids <- c(
  "CRR456609", "CRR456610", "CRR456611", "CRR456612", "CRR456613", "CRR456614",
  "CRR456615", "CRR456616", "CRR456617", "CRR456618", "CRR456619", "CRR456620",
  "CRR456621", "CRR456622", "CRR456623", "CRR456624", "CRR456625", "CRR456626",
  "CRR456627", "CRR456628", "CRR456629", "CRR456630", "CRR456631", "CRR456632"
)

# -----------------------------------------------------------------------------
# 1) Species TPM: first column species name, rest samples (bioproject_runid)
#    Extract runid = suffix after last underscore
# -----------------------------------------------------------------------------
cat("读取菌种 TPM...\n")
dt_species <- fread(file_species_tpm, header = TRUE, sep = "\t")
species_names_raw <- dt_species[[1]]
species_names     <- strip_taxon_prefix(species_names_raw)
sample_cols_species <- names(dt_species)[-1]
runid_from_species  <- sub("^.*_", "", sample_cols_species)

# -----------------------------------------------------------------------------
# 2) Phage TPM: first column vOTU id, rest samples (runid_TPM)
#    Extract runid = strip _TPM suffix
# -----------------------------------------------------------------------------
cat("读取噬菌体 TPM...\n")
dt_phage <- fread(file_phage_tpm, header = TRUE, sep = "\t")
votu_ids   <- dt_phage[[1]]
sample_cols_phage <- names(dt_phage)[-1]
runid_from_phage <- sub("_TPM$", "", sample_cols_phage)

# -----------------------------------------------------------------------------
# 3) Match samples: same runid columns, exclude specified runids
# -----------------------------------------------------------------------------
common_runid <- intersect(runid_from_species, runid_from_phage)
common_runid <- common_runid[!common_runid %in% exclude_runids]
cat("Number of matching common runIDs (excluding the specified runID): ", length(common_runid), "\n")
idx_species_col <- match(common_runid, runid_from_species)
idx_phage_col   <- match(common_runid, runid_from_phage)
# Column indices (1-based, incl. first): species name col 1, samples from 2
j_species <- idx_species_col + 1L
j_phage   <- idx_phage_col + 1L

# Numeric matrix (rows=species/vOTU, cols=matched samples)
Mat_species <- as.matrix(dt_species[, j_species, with = FALSE])
mode(Mat_species) <- "numeric"
rownames(Mat_species) <- species_names

Mat_phage <- as.matrix(dt_phage[, j_phage, with = FALSE])
mode(Mat_phage) <- "numeric"
rownames(Mat_phage) <- votu_ids

n_samples <- length(common_runid)
cat("Sample ", n_samples, "\n")

# -----------------------------------------------------------------------------
# 4) Virulent vOTU set: dual-tool consistent
# -----------------------------------------------------------------------------
bacphlip <- fread(file_bacphlip, header = TRUE, sep = "\t")
if (names(bacphlip)[1] == "" || names(bacphlip)[1] == "V1") names(bacphlip)[1] <- "vOTU"
phatyp <- fread(file_phatyp, header = TRUE)
common_ids <- intersect(bacphlip$vOTU, phatyp$Contig)
bacphlip_common <- bacphlip[vOTU %in% common_ids]
bacphlip_common[, virulent_b := (Virulent > Temperate & Virulent > 0.95)]
phatyp_common <- phatyp[Contig %in% common_ids]
phatyp_common[, virulent_p := (Score > 0.95 & tolower(Pred) == "virulent")]
merge_common <- merge(
  bacphlip_common[virulent_b == TRUE, .(vOTU)],
  phatyp_common[virulent_p == TRUE, .(vOTU = Contig)],
  by = "vOTU"
)
virulent_ids <- merge_common$vOTU
cat("Number of Virulent vOTUs consistent across both software tools:", length(virulent_ids), "\n")

# -----------------------------------------------------------------------------
# 5) Host–vOTU list (species): which vOTUs infect each species
# -----------------------------------------------------------------------------
host_list <- fread(file_host_list, header = TRUE)
host_list[, species := strip_taxon_prefix(host_taxon)]

# Keep vOTUs present in phage matrix (avoid column mismatch)
votu_in_matrix <- rownames(Mat_phage)
parse_virus_list <- function(s) {
  ids <- trimws(strsplit(s, "; ", fixed = TRUE)[[1]])
  ids[ids %in% votu_in_matrix]
}

# -----------------------------------------------------------------------------
# 6) Per species: Spearman between species TPM and total virulent phage TPM
#    Only species in TPM matrix with ≥1 virulent vOTU
# -----------------------------------------------------------------------------
if (n_samples < 3L) stop("The number of matched samples is less than 3; correlation analysis cannot be performed.")

results <- list()
for (i in seq_len(nrow(host_list))) {
  sp <- host_list$species[i]
  virus_list <- parse_virus_list(host_list$virus_id_list[i])
  virulent_infecting <- intersect(virus_list, virulent_ids)
  if (length(virulent_infecting) == 0) next
  if (!sp %in% rownames(Mat_species)) next

  y_species  <- as.numeric(Mat_species[sp, ])
  rows_phage <- which(rownames(Mat_phage) %in% virulent_infecting)
  y_phage_sum <- colSums(Mat_phage[rows_phage, , drop = FALSE], na.rm = TRUE)

  ok <- is.finite(y_species) & is.finite(y_phage_sum)
  if (sum(ok) < 3L) next
  if (sd(y_species[ok]) < 1e-10 || sd(y_phage_sum[ok]) < 1e-10) next

  ct_sp <- cor.test(y_species[ok], y_phage_sum[ok], method = "spearman", exact = FALSE)
  ct_pe <- cor.test(y_species[ok], y_phage_sum[ok], method = "pearson")
  results[[length(results) + 1L]] <- data.table(
    species = sp,
    n_samples = sum(ok),
    n_virulent_vOTUs = length(virulent_infecting),
    spearman_rho = as.numeric(ct_sp$estimate),
    spearman_p = ct_sp$p.value,
    pearson_r = as.numeric(ct_pe$estimate),
    pearson_p = ct_pe$p.value
  )
}

if (length(results) == 0L) {
  stop("No strain simultaneously met the criteria of being present in the species TPM and being infected by at least one lytic vOTU.")
}

res_dt <- rbindlist(results)

# -----------------------------------------------------------------------------
# 7) Multiple testing: FDR (BH) and Bonferroni (Spearman and Pearson)
# -----------------------------------------------------------------------------
res_dt[, p_adj_BH := p.adjust(spearman_p, method = "BH")]
res_dt[, p_adj_Bonferroni := p.adjust(spearman_p, method = "bonferroni")]
res_dt[, pearson_p_adj_BH := p.adjust(pearson_p, method = "BH")]
res_dt[, pearson_p_adj_Bonferroni := p.adjust(pearson_p, method = "bonferroni")]
setorder(res_dt, spearman_p)

cat("FDR (BH) adjusted p_adj_BH < 0.05: ", sum(res_dt$p_adj_BH < 0.05, na.rm = TRUE), "\n")
cat("Bonferroni adjusted p_adj_Bonferroni < 0.05: ", sum(res_dt$p_adj_Bonferroni < 0.05, na.rm = TRUE), "\n")

# -----------------------------------------------------------------------------
# 8) Output: correlation table (Spearman + Pearson)
# -----------------------------------------------------------------------------
fwrite(res_dt, "species_TPM_vs_Virulent_phage_TPM_correlation.csv")
cat("\n Saved: species_TPM_vs_Virulent_phage_TPM_correlation.csv\n")
cat("species, n_samples, n_virulent_vOTUs, spearman_rho, spearman_p, pearson_r, pearson_p, p_adj_BH, p_adj_Bonferroni, pearson_p_adj_BH, pearson_p_adj_Bonferroni\n")

# -----------------------------------------------------------------------------
# 8b) Significant negative species: per-sample species and phage TPM
# -----------------------------------------------------------------------------
sig_neg_species <- res_dt[p_adj_BH < 0.05 & spearman_rho < 0, species]
if (length(sig_neg_species) > 0L) {
  per_sample_list <- list()
  for (sp in sig_neg_species) {
    i <- which(host_list$species == sp)[1]
    if (length(i) == 0) next
    virus_list <- parse_virus_list(host_list$virus_id_list[i])
    virulent_infecting <- intersect(virus_list, virulent_ids)
    if (length(virulent_infecting) == 0 || !sp %in% rownames(Mat_species)) next
    y_species   <- as.numeric(Mat_species[sp, ])
    rows_phage  <- which(rownames(Mat_phage) %in% virulent_infecting)
    y_phage_sum <- colSums(Mat_phage[rows_phage, , drop = FALSE], na.rm = TRUE)
    per_sample_list[[length(per_sample_list) + 1L]] <- data.table(
      species = sp,
      sample_id = common_runid,
      species_TPM = y_species,
      virulent_phage_TPM = y_phage_sum
    )
  }
  if (length(per_sample_list) > 0L) {
    per_sample_dt <- rbindlist(per_sample_list)
    fwrite(per_sample_dt, "species_TPM_virulent_phage_TPM_per_sample_sig_neg.csv")
  }
} else {
  cat("No bacterial species showed a significant negative correlation after FDR correction.\n")
}

# -----------------------------------------------------------------------------
# 9) Binscatter: S. gallolyticus and C. difficile vs virulent phage abundance (negative correlation)
#    X binned by quantiles; Y = mean ± SE of Species TPM [log10(1+x)] per bin
# -----------------------------------------------------------------------------
n_bins <- 10L

plot_binscatter <- function(sp_name, y_species, y_phage_sum) {
  d <- data.frame(
    virulent_phage_TPM = as.numeric(y_phage_sum),
    species_TPM = as.numeric(y_species)
  )
  d <- d[is.finite(d$virulent_phage_TPM) & is.finite(d$species_TPM), ]
  if (nrow(d) < n_bins) return(NULL)
  d$log10_phage <- log10(1 + d$virulent_phage_TPM)
  d$log10_species <- log10(1 + d$species_TPM)
  # Bin X: quantiles first; equal-width if duplicate breaks
  probs <- seq(0, 1, length.out = n_bins + 1L)
  breaks <- quantile(d$log10_phage, probs = probs)
  if (any(duplicated(breaks))) {
    xr <- range(d$log10_phage)
    breaks <- seq(xr[1], xr[2], length.out = n_bins + 1L)
  }
  breaks[1] <- breaks[1] - 1e-9
  breaks[length(breaks)] <- breaks[length(breaks)] + 1e-9
  d$bin <- cut(d$log10_phage, breaks = breaks, include.lowest = TRUE)
  agg <- aggregate(
    cbind(log10_phage, log10_species) ~ bin,
    data = d,
    FUN = function(x) c(mean = mean(x), se = sd(x) / sqrt(length(x)), n = length(x))
  )
  # aggregate columns: col 1 mean, col 2 se, col 3 n
  agg$x_center <- agg$log10_phage[, 1]
  agg$y_mean   <- agg$log10_species[, 1]
  agg$y_se     <- agg$log10_species[, 2]
  agg$n        <- as.integer(agg$log10_species[, 3])
  ct <- cor.test(d$log10_phage, d$log10_species, method = "spearman", exact = FALSE)
  rho <- round(ct$estimate, 3)
  pval <- ct$p.value
  p_lab <- if (pval < 0.001) "p < 0.001" else paste0("p = ", format.pval(pval, digits = 2))
  lbl <- sprintf("Spearman rho = %s\n%s", rho, p_lab)
  p <- ggplot(agg, aes(x = x_center, y = y_mean)) +
    geom_line(colour = "#2D6A4F", linewidth = 1) +
    geom_point(size = 3.5, colour = "#5B9A6E", fill = "white", stroke = 1, shape = 21) +
    geom_errorbar(aes(ymin = y_mean - y_se, ymax = y_mean + y_se), colour = "#2D6A4F", linewidth = 0.6, width = 0.03) +
    annotate("label", x = Inf, y = Inf, label = lbl, hjust = 1.05, vjust = 1.1, size = 4,
             fill = "white", colour = "black", label.padding = unit(0.4, "lines")) +
    labs(
      x = "Virulent vOTU TPM (sum, infecting host) [log10(1+x)]",
      y = "Species TPM [log10(1+x)] (mean +/- SE per bin)",
      title = sp_name
    ) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  return(p)
}

target_species <- c("Streptococcus gallolyticus", "Clostridioides difficile")
for (sp in target_species) {
  i <- which(host_list$species == sp)[1]
  if (length(i) == 0) next
  virus_list <- parse_virus_list(host_list$virus_id_list[i])
  virulent_infecting <- intersect(virus_list, virulent_ids)
  if (length(virulent_infecting) == 0 || !sp %in% rownames(Mat_species)) next
  y_species   <- Mat_species[sp, ]
  rows_phage  <- which(rownames(Mat_phage) %in% virulent_infecting)
  y_phage_sum <- colSums(Mat_phage[rows_phage, , drop = FALSE], na.rm = TRUE)
  p <- plot_binscatter(sp, y_species, y_phage_sum)
  if (!is.null(p)) {
    fname <- gsub(" ", "_", sp)
    ggsave(paste0("binscatter_species_TPM_vs_virulent_phage_TPM_", fname, ".pdf"), p, width = 5.5, height = 5)
    ggsave(paste0("binscatter_species_TPM_vs_virulent_phage_TPM_", fname, ".png"), p, width = 5.5, height = 5, dpi = 150)
    cat("Saved: binscatter_species_TPM_vs_virulent_phage_TPM_", fname, ".pdf / .png\n", sep = "")
  }
}

# -----------------------------------------------------------------------------
# 10) Forest plot: Spearman rho for two species vs virulent phage
# -----------------------------------------------------------------------------
forest_dt <- res_dt[species %in% target_species]
if (nrow(forest_dt) >= 1L) {
  fwrite(forest_dt, "forest_species_virulent_phage_correlation_data.csv")
  cat("Saved: forest_species_virulent_phage_correlation_data.csv\n")
  forest_dt <- copy(forest_dt)
  forest_dt[, species := factor(species, levels = rev(species))]
  x_lim_right <- max(0.05, min(forest_dt$spearman_rho) * 0.5)
  p_forest <- ggplot(forest_dt, aes(x = spearman_rho, y = species)) +
    geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.5, linetype = "dashed") +
    geom_segment(aes(x = 0, xend = spearman_rho, y = species, yend = species),
                 colour = "#2D6A4F", linewidth = 1.2) +
    geom_point(size = 4, colour = "#2D6A4F", fill = "white", stroke = 1, shape = 21) +
    geom_text(aes(label = sprintf("rho = %s, %s", round(spearman_rho, 3),
                                  ifelse(spearman_p < 0.001, "p < 0.001", paste0("p = ", format.pval(spearman_p, digits = 2)))),
              hjust = -0.05, size = 3.5) +
    scale_x_continuous(limits = c(min(forest_dt$spearman_rho) - 0.05, x_lim_right)) +
    labs(x = "Spearman rho (species TPM vs virulent phage TPM)", y = NULL,
         title = "Correlation with virulent phage abundance") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  pdf("forest_species_virulent_phage_correlation.pdf", width = 6, height = 3)
  print(p_forest)
  dev.off()
  png("forest_species_virulent_phage_correlation.png", width = 6, height = 3, units = "in", res = 150)
  print(p_forest)
  dev.off()
  cat("Saved: forest_species_virulent_phage_correlation.pdf / .png\n")
}

# =============================================================================
# Part 2: recompute from per_sample (signal samples only) + forest + phage-stratified bars
# =============================================================================
file_per_sample <- "species_TPM_virulent_phage_TPM_per_sample_sig_neg.csv"
if (!file.exists(file_per_sample)) {
  cat("\nNo bacterial strains with significant negative correlation, or the per_sample table was not generated.\n")
} else {
  dt <- fread(file_per_sample, header = TRUE)
  if (!all(c("species", "sample_id", "species_TPM", "virulent_phage_TPM") %in% names(dt))) {
    stop("per_sample CSV: species, sample_id, species_TPM, virulent_phage_TPM")
  }
  dt[, species_TPM := as.numeric(species_TPM)]
  dt[, virulent_phage_TPM := as.numeric(virulent_phage_TPM)]
  cat("\nPart 2：read per_sample table: ", nrow(dt), "，Bacterial number: ", uniqueN(dt$species), "\n")
  cat("used sample：species_TPM > 0 or virulent_phage_TPM > 0\n")

  results <- list()
  for (sp in unique(dt$species)) {
    d <- dt[species == sp]
    d <- d[is.finite(species_TPM) & is.finite(virulent_phage_TPM)]
    d <- d[(species_TPM > 0) | (virulent_phage_TPM > 0)]
    if (nrow(d) < 3L) next
    if (sd(d$species_TPM) < 1e-10 || sd(d$virulent_phage_TPM) < 1e-10) next
    ct_sp <- cor.test(d$species_TPM, d$virulent_phage_TPM, method = "spearman", exact = FALSE)
    ct_pe <- cor.test(d$species_TPM, d$virulent_phage_TPM, method = "pearson")
    results[[length(results) + 1L]] <- data.table(
      species = sp,
      n_samples = nrow(d),
      spearman_rho = as.numeric(ct_sp$estimate),
      spearman_p = ct_sp$p.value,
      pearson_r = as.numeric(ct_pe$estimate),
      pearson_p = ct_pe$p.value
    )
  }
  if (length(results) == 0L) {
    cat("No strains met the criteria of having an effective sample size of ≥ 3 and non-zero variance\n")
  } else {
    res_dt <- rbindlist(results)
    res_dt[, p_adj_BH := p.adjust(spearman_p, method = "BH")]
    res_dt[, p_adj_Bonferroni := p.adjust(spearman_p, method = "bonferroni")]
    res_dt[, pearson_p_adj_BH := p.adjust(pearson_p, method = "BH")]
    res_dt[, pearson_p_adj_Bonferroni := p.adjust(pearson_p, method = "bonferroni")]
    setorder(res_dt, spearman_p)
    fwrite(res_dt, "reanalyzed_correlation_from_per_sample_sig_neg.csv")
    cat("Saved: reanalyzed_correlation_from_per_sample_sig_neg.csv\n")

    forest_dt <- copy(res_dt)
    forest_dt[, species := factor(species, levels = rev(species))]
    forest_dt[, label := sprintf("rho = %s, %s", round(spearman_rho, 3),
      ifelse(spearman_p < 0.001, "p < 0.001", paste0("p = ", format.pval(spearman_p, digits = 2))))]
    x_lim_left  <- min(forest_dt$spearman_rho) - 0.05
    x_lim_right <- max(0.05, min(forest_dt$spearman_rho) * 0.5)
    p_forest <- ggplot(forest_dt, aes(x = spearman_rho, y = species)) +
      geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.5, linetype = "dashed") +
      geom_segment(aes(x = 0, xend = spearman_rho, y = species, yend = species),
        colour = "#2D6A4F", linewidth = 0.8) +
      geom_point(size = 2, colour = "#2D6A4F", fill = "white", stroke = 0.6, shape = 21) +
      geom_text(aes(label = label), hjust = -0.02, size = 2.5) +
      scale_x_continuous(limits = c(x_lim_left, x_lim_right)) +
      labs(x = "Spearman rho (species TPM vs virulent phage TPM)", y = NULL,
        title = "Correlation with virulent phage abundance (from per-sample sig. neg. table)") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    plot_height <- min(30, max(4, nrow(forest_dt) * 0.12))
    pdf("forest_species_virulent_phage_correlation_from_per_sample.pdf", width = 7, height = plot_height)
    print(p_forest)
    dev.off()
    png("forest_species_virulent_phage_correlation_from_per_sample.png", width = 7, height = plot_height, units = "in", res = 150)
    print(p_forest)
    dev.off()
    cat("Saved: forest_species_virulent_phage_correlation_from_per_sample.pdf / .png\n")

    target_species_scatter <- c("Acetatifactor intestinalis", "Streptococcus gallolyticus", "Clostridioides difficile")
    tier_labels <- c("Low", "Medium", "High")
    list_sum <- list()
    list_annot <- list()
    for (sp in target_species_scatter) {
      d <- dt[species == sp]
      d <- d[is.finite(species_TPM) & is.finite(virulent_phage_TPM)]
      if (nrow(d) < 15L) next
      q13 <- quantile(d$virulent_phage_TPM, probs = 1/3, na.rm = TRUE)
      q23 <- quantile(d$virulent_phage_TPM, probs = 2/3, na.rm = TRUE)
      d[, phage_tier_label := fcase(
        virulent_phage_TPM <= q13, "Low",
        virulent_phage_TPM <= q23, "Medium",
        default = "High"
      )]
      sum_tier <- d[, .(
        mean_species_TPM = mean(species_TPM, na.rm = TRUE),
        se_species_TPM = sd(species_TPM, na.rm = TRUE) / sqrt(.N),
        n = .N
      ), by = phage_tier_label]
      if (nrow(sum_tier) < 3L) next
      sum_tier[, species := sp]
      sum_tier[, phage_tier_label := factor(phage_tier_label, levels = tier_labels)]
      list_sum[[length(list_sum) + 1L]] <- sum_tier
      rr <- res_dt[species == sp]
      rho_val <- if (nrow(rr) > 0) round(rr$spearman_rho[1], 3) else NA_real_
      p_raw <- if (nrow(rr) > 0) rr$p_adj_BH[1] else NA_real_
      padj_val <- if (is.na(p_raw)) "NA" else if (p_raw < 1e-10) formatC(p_raw, format = "e", digits = 1) else format(round(p_raw, 3), scientific = FALSE)
      list_annot[[length(list_annot) + 1L]] <- data.table(
        species = sp,
        label = sprintf("Spearman rho = %s\np_adj = %s", rho_val, padj_val),
        x_tier = "High"
      )
    }
    if (length(list_sum) > 0L) {
      species_order <- c("Acetatifactor intestinalis", "Streptococcus gallolyticus", "Clostridioides difficile")
      combined_sum <- rbindlist(list_sum)
      combined_sum[, species := factor(species, levels = species_order)]
      combined_annot <- rbindlist(list_annot)
      combined_annot[, species := factor(species, levels = species_order)]
      y_max_by_sp <- combined_sum[, .(y_max = max(mean_species_TPM + se_species_TPM, na.rm = TRUE)), by = species]
      combined_annot <- merge(combined_annot, y_max_by_sp, by = "species")
      combined_annot[, y_pos := y_max * 1.12]
      # Virulent-tier sequential: light to dark coral-red (independent of Fig. 4g)
      fill_cols <- c("Low" = "#FEC5BB", "Medium" = "#E76F51", "High" = "#C1121F")
      p_bars <- ggplot(combined_sum, aes(x = phage_tier_label, y = mean_species_TPM, fill = phage_tier_label)) +
        geom_col(colour = "grey30", linewidth = 0.4) +
        geom_errorbar(aes(ymin = mean_species_TPM - se_species_TPM, ymax = mean_species_TPM + se_species_TPM),
          width = 0.22, linewidth = 0.5, colour = "grey30") +
        scale_fill_manual(values = fill_cols, guide = "none") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
        geom_label(data = combined_annot, aes(x = x_tier, y = y_pos, label = label), inherit.aes = FALSE,
          hjust = 0.5, vjust = 1, size = 2.8, fill = "white", colour = "black", label.size = 0.25,
          label.padding = unit(0.4, "lines"), label.r = unit(0.15, "lines")) +
        facet_wrap(~ species, ncol = 3, scales = "free_y") +
        labs(x = "Virulent phage TPM tier", y = "Mean species TPM (\u00b1 SE)",
          title = "Species TPM by virulent phage TPM tier") +
        coord_cartesian(clip = "off") +
        theme_classic(base_size = 11, base_line_size = 0.5) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12, margin = margin(b = 8)),
          axis.title = element_text(face = "plain", size = 11),
          axis.text = element_text(colour = "black", size = 10),
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          axis.line = element_line(colour = "grey30", linewidth = 0.5),
          strip.background = element_rect(fill = "grey92", colour = "grey70", linewidth = 0.4),
          strip.text = element_text(face = "italic", size = 10.5),
          panel.spacing = unit(1, "lines"),
          plot.margin = margin(6, 10, 6, 6)
        )
      pdf("stratified_by_phage_tier_three_species_from_per_sample.pdf", width = 10, height = 4.5)
      print(p_bars)
      dev.off()
      png("stratified_by_phage_tier_three_species_from_per_sample.png", width = 10, height = 4.5, units = "in", res = 300)
      print(p_bars)
      dev.off()
      cat("Saved: stratified_by_phage_tier_three_species_from_per_sample.pdf / .png\n")
    }
  }
}
cat("\nFinished\n")