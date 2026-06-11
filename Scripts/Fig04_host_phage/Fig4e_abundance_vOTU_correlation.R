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
# Correlation analysis and plots: genus/species abundance vs Virulent/Temperate vOTU counts
# Same framework as abundance_vOTU_correlation.R: mean abundance vs vOTUs infecting each host
# Virulent / Temperate × single tool (BACPHLIP, PhaTYP) and dual-tool consistent (Consistent)
# =============================================================================

library(data.table)
library(ggplot2)

# Working directory
if (dir.exists("D:\\F\\MicrobiomeMeta\\global\\Figure\\Phage")) # setwd: using PHAGE_ROOT above
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path)) tryCatch(setwd(dirname(script_path)), error = function(e) NULL)

# -----------------------------------------------------------------------------
# Input files
# -----------------------------------------------------------------------------
file_genus_tpm    <- "taxonomy.genus.tpm.txt"
file_species_tpm  <- "taxonomy.species.tpm.txt"
if (!file.exists(file_genus_tpm) && file.exists(file.path("Functional", file_genus_tpm))) {
  file_genus_tpm <- file.path("Functional", file_genus_tpm)
}
if (!file.exists(file_species_tpm) && file.exists(file.path("Functional", file_species_tpm))) {
  file_species_tpm <- file.path("Functional", file_species_tpm)
}
file_list_genus   <- "host_infection_vOTU_list_per_host_genus.csv"
file_list_species <- "host_infection_vOTU_list_per_host_species.csv"
file_bacphlip     <- "all_vOTUs_lifestyle_predictions.tsv"
file_phatyp       <- "high_quality_combineFilter.phatyp.csv"

strip_taxon_prefix <- function(x) sub("^[a-z]{1,2}__", "", trimws(as.character(x)))

# -----------------------------------------------------------------------------
# 1) Define six vOTU sets: Virulent/Temperate × BACPHLIP / PhaTYP / Consistent
# -----------------------------------------------------------------------------
bacphlip <- fread(file_bacphlip, header = TRUE, sep = "\t")
if (names(bacphlip)[1] == "" || names(bacphlip)[1] == "V1") names(bacphlip)[1] <- "vOTU"
phatyp <- fread(file_phatyp, header = TRUE)

# Single-tool reliable predictions
virulent_bacphlip  <- bacphlip[Virulent > 0.95, vOTU]
temperate_bacphlip <- bacphlip[Temperate > 0.95, vOTU]
virulent_phatyp    <- phatyp[Score > 0.95 & tolower(Pred) == "virulent",  Contig]
temperate_phatyp  <- phatyp[Score > 0.95 & tolower(Pred) == "temperate", Contig]

# Dual-tool consistent predictions
common_ids <- intersect(bacphlip$vOTU, phatyp$Contig)
bacphlip_common <- bacphlip[vOTU %in% common_ids]
bacphlip_common[, reliable_b := (Virulent > 0.95 | Temperate > 0.95)]
bacphlip_common[, lifestyle_b := ifelse(Virulent > Temperate, "Virulent", "Temperate")]
phatyp_common <- phatyp[Contig %in% common_ids]
phatyp_common[, reliable_p := Score > 0.95]
phatyp_common[, lifestyle_p := paste0(toupper(substring(Pred, 1, 1)), substring(Pred, 2))]
merge_common <- merge(
  bacphlip_common[reliable_b == TRUE, .(vOTU, lifestyle_b)],
  phatyp_common[reliable_p == TRUE, .(vOTU = Contig, lifestyle_p)],
  by = "vOTU"
)
both_agree <- merge_common[lifestyle_b == lifestyle_p]
virulent_consistent  <- both_agree[lifestyle_b == "Virulent",  vOTU]
temperate_consistent <- both_agree[lifestyle_b == "Temperate", vOTU]

# Collect into list: lifestyle_type, source, ids
lifestyle_sets <- list(
  Virulent_BACPHLIP   = virulent_bacphlip,
  Virulent_PhaTYP     = virulent_phatyp,
  Virulent_Consistent = virulent_consistent,
  Temperate_BACPHLIP  = temperate_bacphlip,
  Temperate_PhaTYP    = temperate_phatyp,
  Temperate_Consistent = temperate_consistent
)
for (nm in names(lifestyle_sets)) cat(nm, ": ", length(lifestyle_sets[[nm]]), " vOTUs\n", sep = "")

# -----------------------------------------------------------------------------
# 2) Read abundance (same as abundance_vOTU_correlation.R)
# -----------------------------------------------------------------------------
read_abundance <- function(path) {
  dt <- fread(path, header = TRUE)
  taxon_col <- names(dt)[1]
  j_num <- seq.int(2L, ncol(dt))
  for (j in j_num) set(dt, j = j, value = as.numeric(dt[[j]]))
  dt[, mean_abundance := rowMeans(.SD, na.rm = TRUE), .SDcols = j_num]
  out <- dt[, .(taxon = get(taxon_col), mean_abundance)]
  out[, taxon := trimws(as.character(taxon))]
  out[!is.na(mean_abundance) & is.finite(mean_abundance)]
}

cat("读取丰度文件...\n")
abund_genus   <- read_abundance(file_genus_tpm)
abund_species <- read_abundance(file_species_tpm)
abund_genus[, taxon := strip_taxon_prefix(taxon)]
abund_species[, taxon := strip_taxon_prefix(taxon)]

# -----------------------------------------------------------------------------
# 3) From host lists, count lifestyle vOTUs infecting each host
# -----------------------------------------------------------------------------
count_lifestyle_per_host <- function(path_list, ids) {
  dt <- fread(path_list, header = TRUE)
  dt[, taxon := strip_taxon_prefix(host_taxon)]
  dt[, n_vOTUs := sapply(strsplit(virus_id_list, "; ", fixed = TRUE), function(x) {
    sum(trimws(x) %in% ids, na.rm = TRUE)
  })]
  dt[, .(taxon, n_vOTUs)]
}

# For each level and lifestyle_set, produce (taxon, n_vOTUs)
host_genus   <- list()
host_species <- list()
for (nm in names(lifestyle_sets)) {
  host_genus[[nm]]   <- count_lifestyle_per_host(file_list_genus,   lifestyle_sets[[nm]])
  host_species[[nm]] <- count_lifestyle_per_host(file_list_species, lifestyle_sets[[nm]])
}

# -----------------------------------------------------------------------------
# 4) Merge abundance and vOTU counts; keep taxa with both
# -----------------------------------------------------------------------------
merge_abund_n <- function(abund, host_n) {
  merge(abund, host_n, by = "taxon", all = FALSE)
}

merged_genus   <- list()
merged_species <- list()
for (nm in names(lifestyle_sets)) {
  merged_genus[[nm]]   <- merge_abund_n(abund_genus,   host_genus[[nm]])
  merged_species[[nm]] <- merge_abund_n(abund_species, host_species[[nm]])
}

# -----------------------------------------------------------------------------
# 5) Correlation: Spearman + Pearson (same as abundance_vOTU_correlation.R)
# -----------------------------------------------------------------------------
do_correlation <- function(dat, level_label, lifestyle_label, source_label) {
  dat <- copy(dat)
  dat[, log10_abundance := log10(1 + mean_abundance)]
  n <- nrow(dat)
  if (n < 3L) {
    cat(level_label, " ", lifestyle_label, " ", source_label, ": 样本数不足，跳过\n")
    return(list(stats = NULL, dat = dat))
  }
  sp <- cor.test(dat$n_vOTUs, dat$mean_abundance, method = "spearman", exact = FALSE)
  pe <- cor.test(dat$n_vOTUs, dat$log10_abundance, method = "pearson")
  list(
    stats = data.table(
      level = level_label,
      lifestyle = lifestyle_label,
      source = source_label,
      n_taxa = n,
      spearman_rho = as.numeric(sp$estimate),
      spearman_p = sp$p.value,
      pearson_r = as.numeric(pe$estimate),
      pearson_p = pe$p.value
    ),
    dat = dat
  )
}

cor_stats_list <- list()
for (nm in names(lifestyle_sets)) {
  ls_type <- if (grepl("^Virulent", nm)) "Virulent" else "Temperate"
  src     <- sub("^(Virulent|Temperate)_", "", nm)
  res_g   <- do_correlation(merged_genus[[nm]], "Genus",   ls_type, src)
  res_s   <- do_correlation(merged_species[[nm]], "Species", ls_type, src)
  if (!is.null(res_g$stats)) cor_stats_list[[length(cor_stats_list) + 1L]] <- res_g$stats
  if (!is.null(res_s$stats)) cor_stats_list[[length(cor_stats_list) + 1L]] <- res_s$stats
  merged_genus[[nm]]   <- res_g$dat
  merged_species[[nm]] <- res_s$dat
}
cor_stats <- rbindlist(cor_stats_list)

# -----------------------------------------------------------------------------
# 6) Scatter plots: publication style (same as abundance_vOTU_correlation.R)
# -----------------------------------------------------------------------------
theme_pub_scatter <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_blank(),
      axis.title = element_text(size = base_size, colour = "black"),
      axis.text = element_text(size = base_size - 1, colour = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      axis.ticks = element_line(colour = "black", linewidth = 0.35),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      plot.margin = margin(8, 8, 8, 8)
    )
}

plot_correlation <- function(dat, level_label, lifestyle_label, source_label,
                             y_lab = "Mean abundance (log10(1+TPM))") {
  if (nrow(dat) < 3L) return(NULL)
  dat <- copy(dat)
  dat[, log10_abundance := log10(1 + mean_abundance)]
  dat[, log10_n_vOTUs := log10(1 + n_vOTUs)]
  sp <- cor.test(dat$log10_n_vOTUs, dat$log10_abundance, method = "spearman", exact = FALSE)
  pe <- cor.test(dat$log10_n_vOTUs, dat$log10_abundance, method = "pearson")
  lbl <- sprintf("Spearman rho = %s, p %s\nPearson r = %s, p %s",
    round(sp$estimate, 3), format.pval(sp$p.value, digits = 2),
    round(pe$estimate, 3), format.pval(pe$p.value, digits = 2))
  title <- paste0(level_label, " | ", lifestyle_label, " (", source_label, ")")
  ggplot(dat, aes(x = log10_n_vOTUs, y = log10_abundance)) +
    geom_point(alpha = 0.7, size = 2.5, colour = "#2E86AB", stroke = 0.25) +
    geom_smooth(method = "lm", se = TRUE, colour = "#C73E1D", fill = "#F4A582", alpha = 0.25, linewidth = 0.9) +
    annotate("label", x = Inf, y = Inf, label = lbl, hjust = 1.05, vjust = 1.1, size = 3.2,
             fill = "white", colour = "black", label.size = 0.3, label.padding = unit(0.35, "lines")) +
    scale_x_continuous(
      breaks = log10(c(1, 10, 100, 1000, 10000) + 1),
      labels = c("1", "10", "100", "1000", "10000")
    ) +
    labs(
      x = "Number of vOTUs infecting host (log scale)",
      y = y_lab,
      subtitle = title
    ) +
    theme_pub_scatter() +
    theme(plot.subtitle = element_text(hjust = 0.5, size = 10))
}

# Plot: one figure per (level, lifestyle_set)
for (nm in names(lifestyle_sets)) {
  ls_type <- if (grepl("^Virulent", nm)) "Virulent" else "Temperate"
  p_g <- plot_correlation(merged_genus[[nm]], "Genus", ls_type, sub("^(Virulent|Temperate)_", "", nm))
  p_s <- plot_correlation(merged_species[[nm]], "Species", ls_type, sub("^(Virulent|Temperate)_", "", nm))
  if (!is.null(p_g)) {
    ggsave(paste0("abundance_lifestyle_vOTU_correlation_genus_", nm, ".pdf"), p_g, width = 5.5, height = 5)
    ggsave(paste0("abundance_lifestyle_vOTU_correlation_genus_", nm, ".png"), p_g, width = 5.5, height = 5, dpi = 150)
  }
  if (!is.null(p_s)) {
    ggsave(paste0("abundance_lifestyle_vOTU_correlation_species_", nm, ".pdf"), p_s, width = 5.5, height = 5)
    ggsave(paste0("abundance_lifestyle_vOTU_correlation_species_", nm, ".png"), p_s, width = 5.5, height = 5, dpi = 150)
  }
}

# -----------------------------------------------------------------------------
# 7) Output: correlation table + merged tables per combination (optional)
# -----------------------------------------------------------------------------
fwrite(cor_stats, "abundance_lifestyle_vOTU_correlation_stats.csv")
for (nm in names(lifestyle_sets)) {
  merged_genus[[nm]][, log10_abundance := log10(1 + mean_abundance)]
  merged_species[[nm]][, log10_abundance := log10(1 + mean_abundance)]
  fwrite(merged_genus[[nm]][, .(taxon, mean_abundance, log10_abundance, n_vOTUs)],
         paste0("abundance_lifestyle_vOTU_merged_genus_", nm, ".csv"))
  fwrite(merged_species[[nm]][, .(taxon, mean_abundance, log10_abundance, n_vOTUs)],
         paste0("abundance_lifestyle_vOTU_merged_species_", nm, ".csv"))
}
