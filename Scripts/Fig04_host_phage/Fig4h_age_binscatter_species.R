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
# Feces-only: binscatter for two species; overlay four developmental stages (color by Stage)
# - One facet per species; one binned line per Stage (default 3 bins)
# - Spearman rho and p at line end (raw TPM)
#
# Input:
# - species_TPM_virulent_phage_TPM_per_sample_sig_neg.csv
# - meta_qc/meta_augmented_combined.tsv
#
# Output:
# - binscatter_feces_two_species_overlay_stage.pdf / .png
# =============================================================================

library(data.table)
library(ggplot2)

if (dir.exists("D:\\F\\MicrobiomeMeta\\global\\Figure\\Phage")) # setwd: using PHAGE_ROOT above
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path)) tryCatch(setwd(dirname(script_path)), error = function(e) NULL)

file_per_sample <- "species_TPM_virulent_phage_TPM_per_sample_sig_neg.csv"
file_meta <- file.path("meta_qc", "meta_augmented_combined.tsv")
if (!file.exists(file_per_sample)) stop("No input: ", file_per_sample)
if (!file.exists(file_meta)) stop("No input: ", file_meta)

dt <- fread(file_per_sample)
need_cols <- c("species", "sample_id", "species_TPM", "virulent_phage_TPM")
if (!all(need_cols %in% names(dt))) stop("per_sample should comtain column: ", paste(need_cols, collapse = ", "))
dt[, species_TPM := as.numeric(species_TPM)]
dt[, virulent_phage_TPM := as.numeric(virulent_phage_TPM)]

meta <- fread(file_meta)
meta_need <- c("sample_id", "Age_group", "Gut_harmonized", "Gut_location_y", "Gut_location_x")
if (!all(meta_need %in% names(meta))) stop("meta shoudl contain column: ", paste(meta_need, collapse = ", "))
meta <- meta[, ..meta_need]

is_feces <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x %in% c("feces", "faeces")
}
meta[, is_feces := FALSE]
meta[is_feces(Gut_harmonized) | is_feces(Gut_location_y) | is_feces(Gut_location_x), is_feces := TRUE]

dtm <- merge(dt, meta[, .(sample_id, Age_group, is_feces)], by = "sample_id", all.x = TRUE)
dtm <- dtm[is_feces == TRUE]

# Age_group -> 4 stages
dtm[, Stage := fifelse(
  Age_group == "Lactation (0-21d)", "Lactation",
  fifelse(Age_group == "Nursery (22-63d)", "Weaning",
    fifelse(Age_group == "Growing (64-119d)", "Growing",
      fifelse(Age_group %chin% c("Early finishing (120-180d)", "Late finishing (>180d)"), "Adult", NA_character_)
    )
  )
)]
dtm <- dtm[!is.na(Stage)]
dtm[, Stage := factor(Stage, levels = c("Lactation", "Weaning", "Growing", "Adult"))]

target_species <- c("Acetatifactor intestinalis", "CAG-103 sp900317855")
dtm <- dtm[species %in% target_species]
dtm[, species := factor(species, levels = target_species)]
if (nrow(dtm) == 0L) stop("No record")

dtm[, x := log10(1 + species_TPM)]
dtm[, y := log10(1 + virulent_phage_TPM)]

n_bins <- 3L

bin_one_group <- function(xv, yv, n_bins = 3L) {
  ok <- is.finite(xv) & is.finite(yv)
  xv <- xv[ok]; yv <- yv[ok]
  n <- length(xv)
  if (n < max(15L, n_bins * 5L)) return(NULL)
  probs <- seq(0, 1, length.out = n_bins + 1L)
  br <- as.numeric(quantile(xv, probs = probs, na.rm = TRUE))
  if (any(duplicated(br))) {
    xr <- range(xv, na.rm = TRUE)
    if (!is.finite(xr[1]) || !is.finite(xr[2]) || xr[1] == xr[2]) return(NULL)
    br <- seq(xr[1], xr[2], length.out = n_bins + 1L)
  }
  br[1] <- br[1] - 1e-9
  br[length(br)] <- br[length(br)] + 1e-9
  bin <- cut(xv, breaks = br, include.lowest = TRUE)
  dtb <- data.table(x = xv, y = yv, bin = bin)
  dtb[, .(
    x_center = mean(x, na.rm = TRUE),
    y_mean = mean(y, na.rm = TRUE),
    y_se = sd(y, na.rm = TRUE) / sqrt(.N),
    n_bin = .N
  ), by = bin][order(x_center)]
}

spearman_one <- function(xv, yv) {
  ok <- is.finite(xv) & is.finite(yv)
  xv <- xv[ok]; yv <- yv[ok]
  if (length(xv) < 5L) return(list(rho = NA_real_, p = NA_real_, n = length(xv)))
  if (sd(xv) < 1e-10 || sd(yv) < 1e-10) return(list(rho = NA_real_, p = NA_real_, n = length(xv)))
  ct <- suppressWarnings(cor.test(xv, yv, method = "spearman", exact = FALSE))
  list(rho = as.numeric(ct$estimate), p = ct$p.value, n = length(xv))
}

fmt_p <- function(p) {
  p_out <- rep(NA_character_, length(p))
  p_out[is.na(p)] <- "p = NA"
  idx <- !is.na(p) & p < 0.001
  p_out[idx] <- "p < 0.001"
  idx2 <- !is.na(p) & p >= 0.001
  p_out[idx2] <- paste0("p = ", format.pval(p[idx2], digits = 2))
  p_out
}

# Binned data
bins_list <- list()
for (sp in levels(dtm$species)) {
  for (st in levels(dtm$Stage)) {
    dd <- dtm[species == sp & Stage == st]
    b <- bin_one_group(dd$x, dd$y, n_bins = n_bins)
    if (is.null(b)) next
    b[, `:=`(species = sp, Stage = st)]
    bins_list[[length(bins_list) + 1L]] <- b
  }
}
bins <- rbindlist(bins_list, fill = TRUE)
bins[, species := factor(species, levels = target_species)]
bins[, Stage := factor(Stage, levels = levels(dtm$Stage))]

# Label rho/p/n at rightmost point of each line
stats <- dtm[, {
  s <- spearman_one(species_TPM, virulent_phage_TPM)
  .(rho = s$rho, p = s$p, n = s$n)
}, by = .(species, Stage)]
stats[, label := paste0(
  "rho=", ifelse(is.na(rho), "NA", sprintf("%.3f", rho)),
  ", ", fmt_p(p),
  ", n=", n
)]
label_pos <- bins[, .SD[.N], by = .(species, Stage)][, .(species, Stage, x_pos = x_center, y_pos = y_mean)]
stats <- merge(stats, label_pos, by = c("species", "Stage"), all.x = TRUE)
stats[, x_pos := x_pos + 0.06]  # Shift right slightly to avoid overlapping line

stage_cols <- c(
  "Lactation" = "#2b8cbe",
  "Weaning" = "#a6bddb",
  "Growing" = "#fd8d3c",
  "Adult" = "#e6550d"
)

p <- ggplot(bins, aes(x = x_center, y = y_mean, colour = Stage, group = Stage)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.6, fill = "white", stroke = 0.7, shape = 21) +
  geom_errorbar(aes(ymin = y_mean - y_se, ymax = y_mean + y_se), width = 0.03, linewidth = 0.5) +
  geom_text(data = stats[!is.na(x_pos) & !is.na(y_pos)],
    aes(x = x_pos, y = y_pos, label = label, colour = Stage),
    inherit.aes = FALSE, hjust = 0, vjust = 0.5, size = 3) +
  scale_colour_manual(values = stage_cols) +
  facet_wrap(~ species, ncol = 1, scales = "free_y") +
  labs(
    x = "Species TPM [log10(1+x)] (binned)",
    y = "Virulent phage TPM [log10(1+x)] (mean ± SE per bin)",
    colour = "Stage",
    title = "Feces-only: binned trends across stages (overlay; 3 bins)",
    subtitle = "Stage definitions: Lactation 0–21d; Weaning 22–63d; Growing 64–119d; Adult ≥120d"
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12, margin = margin(b = 8)),
    plot.subtitle = element_text(hjust = 0.5, size = 10, margin = margin(b = 6)),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10, colour = "black"),
    axis.line = element_line(colour = "grey30", linewidth = 0.5),
    legend.position = "top",
    strip.background = element_rect(fill = "grey92", colour = "grey70", linewidth = 0.4),
    strip.text = element_text(face = "italic", size = 10.5),
    panel.spacing = unit(1, "lines"),
    plot.margin = margin(6, 28, 6, 6)
  )

ggsave("binscatter_feces_two_species_overlay_stage.pdf", p, width = 8.8, height = 7.2)
ggsave("binscatter_feces_two_species_overlay_stage.png", p, width = 8.8, height = 7.2, dpi = 300)

cat("Saved: binscatter_feces_two_species_overlay_stage.pdf / .png\n")
cat("finished\n")

