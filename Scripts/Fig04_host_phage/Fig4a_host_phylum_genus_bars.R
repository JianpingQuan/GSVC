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
# Genome host prediction: phylum/genus/species vOTU distribution, Top 25 bar plots
# All vOTUs from high_quality_viral_id.txt
# Host unknown = part 1 (not in predictions) + part 2 (in genome but no valid annotation at that level)
# Host annotation normalization: e.g. merge Prevotella copri_A / copri_B / copri_E to Prevotella copri (strip trailing _A/_B/_C etc.)
# Y-axis = host taxonomy, X-axis = log10(vOTU count)
# =============================================================================

library(data.table)
library(ggplot2)

# Working directory (adjust if needed)
if (dir.exists("D:\\F\\MicrobiomeMeta\\global\\Figure\\Phage")) # setwd: using PHAGE_ROOT above
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path)) tryCatch(setwd(dirname(script_path)), error = function(e) NULL)

file_genome   <- "high_quality_Host_prediction_to_genome_m90.csv"
file_all_vOTU <- "high_quality_viral_id.txt"

# Read full high-quality vOTU list
all_vOTU <- fread(file_all_vOTU, header = FALSE, fill = TRUE)[[1]]
all_vOTU <- unique(trimws(all_vOTU))
cat("All high-quality vOTUs with unique quantities (after deduplication):", length(all_vOTU), "\n")

# Read genome host prediction file
d <- fread(file_genome, header = TRUE)
names(d)[1] <- "virus_id"
names(d)[3] <- "host_taxonomy"

# Parse phylum, genus, species from taxonomy string (strip p__/g__/s__ prefixes)
parse_tax <- function(tax) {
  if (is.na(tax) || tax == "") return(list(phylum = NA_character_, genus = NA_character_, species = NA_character_))
  x <- strsplit(tax, ";", fixed = TRUE)[[1]]
  get_level <- function(prefix) {
    idx <- which(grepl(paste0("^", prefix), x))
    if (length(idx) == 0) return(NA_character_)
    sub(paste0("^", prefix), "", x[max(idx)])
  }
  list(
    phylum  = get_level("p__"),
    genus   = get_level("g__"),
    species = get_level("s__")
  )
}

tmp <- lapply(d$host_taxonomy, parse_tax)
d[, phylum  := sapply(tmp, `[[`, "phylum")]
d[, genus   := sapply(tmp, `[[`, "genus")]
d[, species := sapply(tmp, `[[`, "species")]

# Drop empty or NA species (retain phylum, genus)
d[species == "" | is.na(species), species := NA_character_]

# Host annotation normalization: merge Prevotella copri_A / copri_B / copri_E etc. to Prevotella copri (strip trailing _A/_B/_C etc.)
# Use merged names consistently in plots and statistics
normalize_taxon <- function(x) {
  if (is.na(x) || !nzchar(trimws(x))) return(x)
  sub("_[A-Z][0-9]*$", "", x)  # Strip trailing _A, _B, _E, _A1, _B2 etc.
}
d[, phylum  := sapply(phylum, normalize_taxon)]
d[, genus   := sapply(genus,  normalize_taxon)]
d[, species := sapply(species, normalize_taxon)]
# Normalization may yield empty strings; clean again
d[phylum == "" | is.na(phylum), phylum := NA_character_]
d[genus == "" | is.na(genus), genus := NA_character_]
d[species == "" | is.na(species), species := NA_character_]

# Unify phylum name: Bacillota → Firmicutes (common in high-impact journals)
d[phylum == "Bacillota", phylum := "Firmicutes"]

# Host unknown has two parts:
# Part 1: vOTUs not in genome predictions = all high-quality viral ids − viral ids with host in genome
# Part 2: vOTUs in genome file but without valid annotation at that taxonomic level
genome_virus_ids <- unique(d$virus_id)
not_in_genome    <- setdiff(all_vOTU, genome_virus_ids)
n_part1          <- length(not_in_genome)   # vOTUs not in prediction results
cat("Predicted number of host vOTUs (in the genome file):", length(genome_virus_ids), "\n")
cat("Number of vOTUs not in the predicted results (Host unknown, Part 1):", n_part1, "\n")

# Count by phylum/genus/species: annotated vOTUs by taxonomy; Host unknown = part 1 + part 2
count_vOTU_with_unknown <- function(dt, level_col, all_ids, genome_ids, n_not_in_genome) {
  u <- unique(dt[!is.na(get(level_col)) & nzchar(trimws(get(level_col))), .(virus_id, taxon = get(level_col))])
  by_taxon <- u[, .(n_vOTU = .N), by = taxon]
  with_annot <- unique(u$virus_id)
  # Part 2: in genome but no valid annotation at this level
  in_genome_no_annot <- setdiff(genome_ids, with_annot)
  n_part2 <- length(in_genome_no_annot)
  n_unknown <- n_not_in_genome + n_part2
  rbind(by_taxon, data.table(taxon = "Host unknown", n_vOTU = n_unknown))
}

cat("Tabulate vOTU counts at the phylum, genus, and species levels (including "Host unknown" = Part 1 + Part 2)...\n")
phylum_n  <- count_vOTU_with_unknown(d, "phylum",  all_vOTU, genome_virus_ids, n_part1)
genus_n   <- count_vOTU_with_unknown(d, "genus",   all_vOTU, genome_virus_ids, n_part1)
species_n <- count_vOTU_with_unknown(d, "species", all_vOTU, genome_virus_ids, n_part1)

# Host unknown counts per level and two-part breakdown (part 2 differs by level)
host_unknown_breakdown <- data.table(
  level   = c("phylum", "genus", "species"),
  part1_not_in_genome = rep(n_part1, 3),
  part2_in_genome_no_annot = c(
    phylum_n[taxon == "Host unknown", n_vOTU] - n_part1,
    genus_n[taxon == "Host unknown", n_vOTU] - n_part1,
    species_n[taxon == "Host unknown", n_vOTU] - n_part1
  )
)
host_unknown_breakdown[, host_unknown_total := part1_not_in_genome + part2_in_genome_no_annot]
for (i in 1:3) {
  cat(host_unknown_breakdown[i, level], "level Host unknown:", host_unknown_breakdown[i, host_unknown_total],
      " (Part 1: Not in the predicted results =", host_unknown_breakdown[i, part1_not_in_genome],
      "; Part 2: Present in the genome but lacking annotation at that level =", host_unknown_breakdown[i, part2_in_genome_no_annot], ")\n")
}

# Top 25 predicted taxa + Host unknown (up to 26 items)
top_n <- 25
phylum_top  <- rbind(phylum_n[taxon != "Host unknown"][order(-n_vOTU)][1:top_n], phylum_n[taxon == "Host unknown"])
phylum_top  <- phylum_top[order(-n_vOTU)]
genus_top   <- rbind(genus_n[taxon != "Host unknown"][order(-n_vOTU)][1:top_n], genus_n[taxon == "Host unknown"])
genus_top   <- genus_top[order(-n_vOTU)]
species_top <- rbind(species_n[taxon != "Host unknown"][order(-n_vOTU)][1:top_n], species_n[taxon == "Host unknown"])
species_top <- species_top[order(-n_vOTU)]

# Phylum for genus/species (for coloring): most frequent phylum per taxon
phylum_by_genus   <- d[!is.na(genus) & nzchar(trimws(genus)), .N, by = .(genus, phylum)][, .(phylum = phylum[which.max(N)]), by = genus]
phylum_by_species <- d[!is.na(species) & nzchar(trimws(species)), .N, by = .(species, phylum)][, .(phylum = phylum[which.max(N)]), by = species]
setnames(phylum_by_genus, "genus", "taxon")
setnames(phylum_by_species, "species", "taxon")
genus_top   <- merge(genus_top, phylum_by_genus, by = "taxon", all.x = TRUE)
genus_top[is.na(phylum) | taxon == "Host unknown", phylum := "Host unknown"]
species_top <- merge(species_top, phylum_by_species, by = "taxon", all.x = TRUE)
species_top[is.na(phylum) | taxon == "Host unknown", phylum := "Host unknown"]

# Publication theme (high-impact journal style): white background, clear axes, readable fonts, minimal grid
theme_pub <- function(base_size = 11, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title    = element_blank(),
      axis.title    = element_text(size = base_size, colour = "black"),
      axis.text     = element_text(size = base_size - 1, colour = "black"),
      axis.line     = element_line(colour = "black", linewidth = 0.35),
      axis.ticks    = element_line(colour = "black", linewidth = 0.35),
      legend.title  = element_text(size = base_size, face = "bold"),
      legend.text   = element_text(size = base_size - 1),
      legend.key.size = unit(4, "mm"),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      plot.margin    = margin(6, 8, 6, 6)
    )
}

# AAAS-style phylum colors (match panel b); panel a uses red→blue gradient by vOTU rank
AAAS_COLORS <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F")

# Phylum level: gradient (red→white→blue, high count=red), Host unknown gray, no legend
bar_plot_phylum <- function(dat, xlab = "log10(vOTU count)") {
  dat <- copy(dat)
  setorder(dat, -n_vOTU)   # Sort descending by count; largest bar at top when plotted
  dat[, taxon := factor(taxon, levels = rev(dat$taxon))]  # After flip: first level at bottom, last at top
  dat[, log_n := log10(pmax(n_vOTU, 1))]
  dat[, is_unknown := (taxon == "Host unknown")]
  n_known <- sum(!dat$is_unknown)
  grad_cols <- colorRampPalette(c(AAAS_COLORS[1], "#F7F7F7", AAAS_COLORS[4]))(max(n_known, 1L))
  dat[, fill_col := "gray70"]
  dat[is_unknown == FALSE, fill_col := grad_cols[seq_len(n_known)]]
  ggplot(dat, aes(x = taxon, y = log_n, fill = fill_col)) +
    geom_col(width = 0.72, linewidth = 0.2, colour = "grey92") +
    scale_fill_identity(guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
    coord_flip() +
    labs(y = xlab, x = "Host taxonomy") +
    theme_pub()
}
p_phylum <- bar_plot_phylum(phylum_top)

# Genus/species level: color by phylum with legend; bars sorted top-to-bottom by count (largest on top)
bar_plot_by_phylum <- function(dat, xlab = "log10(vOTU count)") {
  dat <- copy(dat)
  setorder(dat, -n_vOTU)
  dat[, taxon := factor(taxon, levels = rev(dat$taxon))]
  dat[, log_n := log10(pmax(n_vOTU, 1))]
  phyla_ord <- setdiff(unique(dat$phylum), "Host unknown")
  n_phyla <- length(phyla_ord)
  pal_phyla <- if (n_phyla > 0) setNames(
    colorRampPalette(AAAS_COLORS)(n_phyla),
    phyla_ord
  ) else character(0)
  pal <- c(pal_phyla, "Host unknown" = "gray70")
  ggplot(dat, aes(x = taxon, y = log_n, fill = phylum)) +
    geom_col(width = 0.72, linewidth = 0.2, colour = "grey92") +
    scale_fill_manual(values = pal, name = "Phylum", na.value = "gray70") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
    coord_flip() +
    labs(y = xlab, x = "Host taxonomy") +
    theme_pub() +
    theme(legend.position = "right", axis.text.y = element_text(face = "italic"))
}
p_genus   <- bar_plot_by_phylum(genus_top)
p_species <- bar_plot_by_phylum(species_top)

# Save full stats (incl. Host unknown) and Top 25 + Host unknown plot data
fwrite(host_unknown_breakdown, "host_unknown_two_parts_breakdown.csv")
fwrite(phylum_n,   "host_taxonomy_vOTU_full_phylum.csv")
fwrite(genus_n,    "host_taxonomy_vOTU_full_genus.csv")
fwrite(species_n,  "host_taxonomy_vOTU_full_species.csv")
fwrite(phylum_top, "host_taxonomy_vOTU_top25_phylum.csv")
fwrite(genus_top,  "host_taxonomy_vOTU_top25_genus.csv")
fwrite(species_top,"host_taxonomy_vOTU_top25_species.csv")

# Save figures: individual panels + three-panel combined
ggsave("host_taxonomy_vOTU_top25_phylum.pdf",  p_phylum,  width = 4, height = 7)
ggsave("host_taxonomy_vOTU_top25_genus.pdf",   p_genus,   width = 6, height = 7)
ggsave("host_taxonomy_vOTU_top25_species.pdf",  p_species, width = 9, height = 7)
if (requireNamespace("gridExtra", quietly = TRUE)) {
  pdf("host_taxonomy_vOTU_top25_combined.pdf", width = 9, height = 16)
  gridExtra::grid.arrange(p_phylum, p_genus, p_species, ncol = 1)
  dev.off()
} else {
  cat("After installing gridExtra, you can save a version with the three plots combined: install.packages('gridExtra')\n")
}
