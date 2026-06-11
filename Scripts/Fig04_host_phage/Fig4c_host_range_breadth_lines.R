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
# Host range and infection load statistics
# 1) Phage host count: distribution of vOTUs infecting 1/2/3/... hosts (phylum/genus/species) → bar plot
# 2) Host phage load: distribution of hosts infected by 1/2/3/... vOTUs → bar plot + table
# All three levels (phylum, genus, species); plot style matches Host range breadth (gradient colors, count and pct labels)
# =============================================================================

library(data.table)
library(ggplot2)

# Working directory
if (dir.exists("D:\\F\\MicrobiomeMeta\\global\\Figure\\Phage")) # setwd: using PHAGE_ROOT above
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(script_path)) tryCatch(setwd(dirname(script_path)), error = function(e) NULL)

file_genome <- "high_quality_Host_prediction_to_genome_m90.csv"

# Read and parse taxonomy
d <- fread(file_genome, header = TRUE)
names(d)[1] <- "virus_id"
names(d)[3] <- "host_taxonomy"
if ("Confidence score" %in% names(d)) d[, confidence := as.numeric(`Confidence score`)]
if ("Main method" %in% names(d)) d[, main_method := as.character(`Main method`)]

parse_tax <- function(tax) {
  if (is.na(tax) || tax == "") return(list(phylum = NA_character_, genus = NA_character_, species = NA_character_))
  x <- strsplit(tax, ";", fixed = TRUE)[[1]]
  get_level <- function(prefix) {
    idx <- which(grepl(paste0("^", prefix), x))
    if (length(idx) == 0) return(NA_character_)
    sub(paste0("^", prefix), "", x[max(idx)])
  }
  list(phylum = get_level("p__"), genus = get_level("g__"), species = get_level("s__"))
}
tmp <- lapply(d$host_taxonomy, parse_tax)
d[, phylum  := sapply(tmp, `[[`, "phylum")]
d[, genus   := sapply(tmp, `[[`, "genus")]
d[, species := sapply(tmp, `[[`, "species")]
d[species == "" | is.na(species), species := NA_character_]

normalize_taxon <- function(x) {
  if (is.na(x) || !nzchar(trimws(x))) return(x)
  sub("_[A-Z][0-9]*$", "", x)
}
d[, phylum  := sapply(phylum, normalize_taxon)]
d[, genus   := sapply(genus,  normalize_taxon)]
d[, species := sapply(species, normalize_taxon)]
d[phylum == "" | is.na(phylum), phylum := NA_character_]
d[genus == "" | is.na(genus), genus := NA_character_]
d[species == "" | is.na(species), species := NA_character_]
d[phylum == "Bacillota", phylum := "Firmicutes"]

# Level names (for labels and filenames)
# capped at 20: bin host counts >20 as 20 to avoid overly long bars
level_names <- list(
  phylum  = list(label = "phylum",  xlab = "Number of host phyla per vOTU",   xlab2 = "Number of vOTUs per host"),
  genus   = list(label = "genus",   xlab = "Number of host genera per vOTU ",  xlab2 = "Number of vOTUs per host"),
  species = list(label = "species", xlab = "Number of host species per vOTU", xlab2 = "Number of vOTUs per host")
)
cap_at <- 20L

# Publication-style theme (no title)
theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_blank(),
      axis.title = element_text(size = base_size, colour = "black"),
      axis.text  = element_text(size = base_size - 1, colour = "black"),
      axis.line  = element_line(colour = "black", linewidth = 0.35),
      axis.ticks = element_line(colour = "black", linewidth = 0.35),
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white"),
      plot.margin = margin(6, 8, 6, 6)
    )
}

# Two gradient palettes: phage→host vs host→phage
# (1) Phage host count: blue→purple→orange→red
gradient_colors_breadth <- colorRampPalette(c("#2166AC", "#67A9CF", "#B2182B", "#EF8A62", "#FDDBC7", "#F7F7F7"))(cap_at)
# (2) Host phage count: green→cyan→purple
gradient_colors_infection <- colorRampPalette(c("#006837", "#41AB5D", "#ADDD8E", "#253494", "#4575B4", "#ABD9E9"))(cap_at)
# X-axis: last bin 20 shown as >=20 (unified bar labels below)
# Note: few phyla; most are infected by many vOTUs, so n_vOTUs often 1–9 or >=20; bins 10–19 often zero (expected)

# -----------------------------------------------------------------------------
# (1) Phage host count: distinct hosts per vOTU → distribution bar plot + per-vOTU host detail
# -----------------------------------------------------------------------------
host_range_plots <- list()
host_range_tables <- list()
host_range_per_vOTU_tables <- list()  # Which hosts (phylum/genus/species) and host count per virus_id

for (lev in c("phylum", "genus", "species")) {
  col <- lev
  meta <- level_names[[lev]]
  # Distinct host count per vOTU (rows with annotation at this level only)
  dd <- d[!is.na(get(col)) & nzchar(trimws(get(col)))]
  if (nrow(dd) == 0) next
  # Per vOTU: host count + host list (comma-separated) for multi-host infections
  per_vOTU <- dd[, .(
    n_hosts = uniqueN(get(col)),
    host_taxa = paste(sort(unique(get(col))), collapse = "; ")
  ), by = virus_id]
  setorder(per_vOTU, -n_hosts, virus_id)
  host_range_per_vOTU_tables[[lev]] <- per_vOTU

  n_hosts_per_vOTU <- dd[, .(n_hosts = uniqueN(get(col))), by = virus_id]
  # Distribution: how many vOTUs have 1 host, 2 hosts, …
  dist <- n_hosts_per_vOTU[, .N, by = n_hosts][order(n_hosts)]
  dist[, n_hosts_cap := pmin(n_hosts, cap_at)]
  dist_cap <- dist[, .(n_vOTUs = sum(N)), by = n_hosts_cap][order(n_hosts_cap)]
  setnames(dist_cap, "n_hosts_cap", "n_hosts")
  total <- sum(dist_cap$n_vOTUs)
  dist_cap[, pct := round(100 * n_vOTUs / total, 1)]
  # Phylum: label count and (pct) on bars; genus/species: pct only
  if (lev == "phylum") {
    dist_cap[, label := paste0(format(n_vOTUs, big.mark = ","), " (", pct, "%)")]
  } else {
    dist_cap[, label := paste0(pct, "%")]
  }
  host_range_tables[[lev]] <- copy(dist_cap)

  # Bar plot: x = n_hosts (1..cap_at), y = n_vOTUs, first gradient; narrower bars at phylum level
  dist_cap[, n_hosts := factor(n_hosts, levels = 1:cap_at)]
  fill_vals <- setNames(gradient_colors_breadth[seq_len(cap_at)], as.character(1:cap_at))
  bar_width <- if (lev == "phylum") 0.5 else 0.7
  p <- ggplot(dist_cap, aes(x = n_hosts, y = n_vOTUs, fill = n_hosts)) +
    geom_col(width = bar_width, linewidth = 0.2, colour = "grey92") +
    geom_text(aes(label = label), vjust = -0.3, size = 3, colour = "black") +
    scale_fill_manual(values = fill_vals, guide = "none") +
    scale_x_discrete(labels = c(as.character(1:(cap_at - 1L)), ">=20")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    coord_cartesian(clip = "off") +
    labs(x = meta$xlab, y = "Number of vOTUs") +
    theme_pub()
  host_range_plots[[lev]] <- p
}

# Save figures and tables for (1)
for (lev in names(host_range_plots)) {
  ggsave(paste0("host_range_breadth_", lev, ".pdf"), host_range_plots[[lev]], width = 7, height = 5)
  ggsave(paste0("host_range_breadth_", lev, ".png"), host_range_plots[[lev]], width = 7, height = 5, dpi = 150)
}
for (lev in names(host_range_tables)) {
  fwrite(host_range_tables[[lev]], paste0("host_range_breadth_", lev, "_table.csv"))
}
for (lev in names(host_range_per_vOTU_tables)) {
  fwrite(host_range_per_vOTU_tables[[lev]], paste0("host_range_per_vOTU_", lev, ".csv"))
}

# -----------------------------------------------------------------------------
# (2) Host phage load: vOTUs per host → bar plot + per-host vOTU count and list
# -----------------------------------------------------------------------------
host_infection_plots <- list()
host_infection_dist_tables <- list()  # Distribution: host count for k vOTUs
host_infection_detail_tables <- list() # Each host and its vOTU count
host_infection_vOTU_list_tables <- list() # virus_id list per host

for (lev in c("phylum", "genus", "species")) {
  col <- lev
  meta <- level_names[[lev]]
  dd <- d[!is.na(get(col)) & nzchar(trimws(get(col)))]
  if (nrow(dd) == 0) next
  # Per host (taxon): vOTU count + infecting virus_id list (semicolon-separated)
  per_host <- dd[, .(
    n_vOTUs = uniqueN(virus_id),
    virus_id_list = paste(sort(unique(virus_id)), collapse = "; ")
  ), by = c(col)]
  setnames(per_host, col, "host_taxon")
  setorder(per_host, -n_vOTUs, host_taxon)
  host_infection_vOTU_list_tables[[lev]] <- per_host
  host_infection_detail_tables[[lev]] <- per_host[, .(host_taxon, n_vOTUs)]

  n_vOTUs_per_host <- dd[, .(n_vOTUs = uniqueN(virus_id)), by = c(col)]
  setnames(n_vOTUs_per_host, col, "host_taxon")

  # Distribution: hosts with 1 vOTU, 2 vOTUs, …
  dist <- n_vOTUs_per_host[, .N, by = n_vOTUs][order(n_vOTUs)]
  dist[, n_vOTUs_cap := pmin(n_vOTUs, cap_at)]
  dist_cap <- dist[, .(n_hosts = sum(N)), by = n_vOTUs_cap][order(n_vOTUs_cap)]
  setnames(dist_cap, "n_vOTUs_cap", "n_vOTUs")
  total <- sum(dist_cap$n_hosts)
  dist_cap[, pct := round(100 * n_hosts / total, 1)]
  dist_cap[, label := paste0(pct, "%")]
  host_infection_dist_tables[[lev]] <- copy(dist_cap)

  # Bar plot: second gradient for host→phage; color bins 1..20 (phylum may have sparse bins)
  dist_cap[, n_vOTUs := factor(n_vOTUs, levels = 1:cap_at)]
  fill_vals <- setNames(gradient_colors_infection[seq_len(cap_at)], as.character(1:cap_at))
  p <- ggplot(dist_cap, aes(x = n_vOTUs, y = n_hosts, fill = n_vOTUs)) +
    geom_col(width = 0.7, linewidth = 0.2, colour = "grey92") +
    geom_text(aes(label = label), vjust = -0.3, size = 3, colour = "black") +
    scale_fill_manual(values = fill_vals, guide = "none") +
    scale_x_discrete(labels = c(as.character(1:(cap_at - 1L)), ">=20")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    coord_cartesian(clip = "off") +
    labs(x = meta$xlab2, y = "Number of hosts") +
    theme_pub()
  host_infection_plots[[lev]] <- p
}

# Save figures and tables for (2)
for (lev in names(host_infection_plots)) {
  ggsave(paste0("host_infection_load_", lev, ".pdf"), host_infection_plots[[lev]], width = 7, height = 5)
  ggsave(paste0("host_infection_load_", lev, ".png"), host_infection_plots[[lev]], width = 7, height = 5, dpi = 150)
}
for (lev in names(host_infection_dist_tables)) {
  fwrite(host_infection_dist_tables[[lev]], paste0("host_infection_load_distribution_", lev, "_table.csv"))
}
for (lev in names(host_infection_detail_tables)) {
  fwrite(host_infection_detail_tables[[lev]], paste0("host_infection_load_per_host_", lev, "_table.csv"))
}
for (lev in names(host_infection_vOTU_list_tables)) {
  fwrite(host_infection_vOTU_list_tables[[lev]], paste0("host_infection_vOTU_list_per_host_", lev, ".csv"))
}

# -----------------------------------------------------------------------------
# (4) Specific hosts (species): species infected by exactly one vOTU
# - Output detail table + summary plot (count unique-infection species by host phylum)
# -----------------------------------------------------------------------------
if (!is.null(host_infection_vOTU_list_tables[["species"]])) {
  dd_sp <- d[!is.na(species) & nzchar(trimws(species))]
  # Per species: unique vOTU count + vOTU list + phylum/genus for aggregation
  per_sp <- dd_sp[, .(
    n_vOTUs = uniqueN(virus_id),
    virus_id_list = paste(sort(unique(virus_id)), collapse = "; "),
    phylum = na.omit(unique(phylum))[1],
    genus = na.omit(unique(genus))[1]
  ), by = .(host_taxon = species)]
  per_sp[is.na(phylum) | !nzchar(trimws(phylum)), phylum := "Unclassified"]
  per_sp[is.na(genus) | !nzchar(trimws(genus)), genus := "Unclassified"]
  setorder(per_sp, host_taxon)
  fwrite(per_sp, "host_species_infection_summary.tsv", sep = "\t")

  uniq_sp <- per_sp[n_vOTUs == 1L]

  # Add unique vOTU confidence/method to uniq_sp (from raw host prediction table d)
  # - One vOTU per species by definition; raw table may have multiple rows; take max confidence
  d_sp <- d[!is.na(species) & nzchar(trimws(species)) & !is.na(confidence)]
  d_sp[, host_taxon := species]
  conf_per_pair <- d_sp[, .(
    confidence_max = max(confidence, na.rm = TRUE),
    confidence_median = as.numeric(stats::median(confidence, na.rm = TRUE)),
    main_method_top = names(sort(table(main_method), decreasing = TRUE))[1]
  ), by = .(host_taxon, virus_id)]
  setnames(conf_per_pair, "virus_id", "virus_id_only")

  uniq_sp[, virus_id_only := tstrsplit(virus_id_list, ";\\s*")[[1]]]
  uniq_sp <- merge(uniq_sp, conf_per_pair, by = c("host_taxon", "virus_id_only"), all.x = TRUE)
  uniq_sp[is.na(confidence_max), confidence_max := NA_real_]
  uniq_sp[is.na(confidence_median), confidence_median := NA_real_]
  uniq_sp[is.na(main_method_top) | !nzchar(trimws(main_method_top)), main_method_top := "NA"]

  # For manual review: unique-infection species list (unique vOTU + confidence)
  fwrite(uniq_sp, "host_species_unique_vOTU_only.tsv", sep = "\t")

  # Plot: species infected by a single vOTU (Top N by confidence_max)
  # - y-axis: Host species (truncated if too long)
  # - x-axis: Confidence score (max)
  # - Point color: Phylum (Top phyla + Other) for readability
  top_n <- 40L
  uniq_plot <- uniq_sp[!is.na(confidence_max)]
  if (nrow(uniq_plot) > 0) {
    uniq_plot <- uniq_plot[order(-confidence_max)]
    if (nrow(uniq_plot) > top_n) {
      uniq_plot <- uniq_plot[1:top_n]
    }
    # Show Top 8 phyla; rest as Other
    ph_ct <- uniq_plot[, .N, by = phylum][order(-N)]
    keep_ph <- ph_ct[1:min(8L, nrow(ph_ct)), phylum]
    uniq_plot[, phylum_plot := ifelse(phylum %in% keep_ph, phylum, "Other")]
    uniq_plot[, host_label := host_taxon]
    uniq_plot[nchar(host_label) > 34, host_label := paste0(substr(host_label, 1, 33), "…")]
    uniq_plot[, host_label := factor(host_label, levels = rev(unique(host_label)))]

    p_species_unique <- ggplot(uniq_plot, aes(x = confidence_max, y = host_label)) +
      geom_segment(aes(x = 0, xend = confidence_max, yend = host_label), colour = "grey85", linewidth = 0.6) +
      geom_point(aes(colour = phylum_plot), size = 2.6, alpha = 0.95) +
      scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), expand = expansion(mult = c(0.01, 0.02))) +
      labs(
        x = "Host assignment confidence score (max)",
        y = NULL,
        colour = "Host phylum",
        subtitle = paste0("Host species infected by exactly one vOTU (Top ", top_n, " by confidence)")
      ) +
      theme_classic(base_size = 11) +
      theme(
        plot.title = element_blank(),
        plot.subtitle = element_text(size = 10, colour = "grey30"),
        axis.text.y = element_text(size = 9),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        legend.position = "right",
        plot.margin = margin(6, 10, 6, 6)
      )

    ggsave("Fig 4x - host_species_unique_vOTU_lollipop_by_confidence.pdf", p_species_unique, width = 8.8, height = 6.2)
    ggsave("Fig 4x - host_species_unique_vOTU_lollipop_by_confidence.png", p_species_unique, width = 8.8, height = 6.2, dpi = 300)
  }

  # Plot: genus-level summary — count of unique-infection species per genus
  # - x-axis: count; y-axis: genus (Top N + Other) for long labels
  uniq_sp_genus <- uniq_sp[, .(n_species_unique = .N), by = genus][order(-n_species_unique)]
  top_g <- 10L
  uniq_sp_genus2 <- uniq_sp_genus[1:min(top_g, nrow(uniq_sp_genus))]
  uniq_sp_genus2[, genus := as.character(genus)]
  uniq_sp_genus2[is.na(genus) | !nzchar(trimws(genus)), genus := "Unclassified"]
  # Sort ascending; with coord_flip, descending top to bottom
  setorder(uniq_sp_genus2, n_species_unique, genus)
  uniq_sp_genus2[, genus := factor(genus, levels = uniq_sp_genus2$genus)]

  # Journal style: horizontal bars + value labels; light grid; Arial; compact single-column size
  maxv <- max(uniq_sp_genus2$n_species_unique, na.rm = TRUE)
  uniq_sp_genus2[, label := format(n_species_unique, big.mark = ",")]

  p_unique_genus_bar <- ggplot(uniq_sp_genus2, aes(x = genus, y = n_species_unique)) +
    geom_col(width = 0.72, fill = "#2F6DAE", colour = "white", linewidth = 0.2) +
    geom_text(
      aes(label = label),
      hjust = -0.15,
      size = 3.4,
      colour = "#1a1a1a"
    ) +
    coord_flip(clip = "off") +
    scale_y_continuous(
      limits = c(0, maxv * 1.12),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      x = NULL,
      y = "No. host species with exactly one vOTU",
      subtitle = paste0("Top ", min(top_g, nrow(uniq_sp_genus)), " genera")
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_text(size = 10, colour = "grey30"),
      axis.title.x = element_text(size = 11),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      axis.ticks = element_line(colour = "black", linewidth = 0.35),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(colour = "#E6E6E6", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      plot.margin = margin(6, 18, 6, 6)
    )

  ggsave("Fig 4x - unique_infected_species_count_by_genus.pdf", p_unique_genus_bar, width = 6.8, height = 4.2)
  ggsave("Fig 4x - unique_infected_species_count_by_genus.png", p_unique_genus_bar, width = 6.8, height = 4.2, dpi = 300)
}

# -----------------------------------------------------------------------------
# (3) Line plots: phylum/genus/species same axes, phage→host and host→phage distributions
# -----------------------------------------------------------------------------
line_theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_blank(),
      axis.title = element_text(size = base_size, colour = "black"),
      axis.text = element_text(size = base_size - 1, colour = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.35),
      axis.ticks = element_line(colour = "black", linewidth = 0.35),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      legend.position = "right",
      legend.key.width = unit(8, "mm"),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      plot.margin = margin(6, 8, 6, 6)
    )
}
level_colors <- c(Phylum = "#2166AC", Genus = "#2CA02C", Species = "#D62728")

# Phage host count distribution: three lines, x = host count (1..>=20), y = pct (%)
rbreadth <- rbindlist(lapply(names(host_range_tables), function(lev) {
  dt <- copy(host_range_tables[[lev]])
  dt[, level := factor(paste0(toupper(substring(lev, 1, 1)), substring(lev, 2)), levels = c("Phylum", "Genus", "Species"))]
  dt[, x := as.integer(n_hosts)]
  dt[, y := n_vOTUs]
  dt[, .(level, x, y)]
}))
rbreadth[, y_pct := 100 * y / sum(y), by = level]
x_breaks <- seq.int(1L, cap_at)
x_breaks_lab <- as.character(x_breaks)
x_breaks_lab[cap_at] <- ">=20"

p_line_breadth <- ggplot(rbreadth, aes(x = x, y = y_pct, colour = level)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_colour_manual(values = level_colors, name = "Level") +
  scale_x_continuous(breaks = x_breaks, labels = x_breaks_lab) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), labels = function(z) paste0(z, "%")) +
  labs(x = "Number of hosts per vOTU", y = "Proportion (%)") +
  line_theme_pub()

# ---- Inset: generalist fraction (hosts >= 2) ----
# Inset with same Level colors to support generalist conclusion in text.
library(grid)
inset_dt <- copy(rbreadth)
inset_dt[, is_generalist := x >= 2L]
inset_sum <- inset_dt[, .(pct_generalist = sum(y[is_generalist]) / sum(y) * 100), by = level]
inset_sum[, pct_specialist := 100 - pct_generalist]
inset_long <- melt(
  inset_sum,
  id.vars = "level",
  measure.vars = c("pct_specialist", "pct_generalist"),
  variable.name = "class",
  value.name = "pct"
)
inset_long[, class := factor(class, levels = c("pct_specialist", "pct_generalist"),
                             labels = c("Specialist (1 host)", "Generalist (>=2 hosts)"))]

p_inset <- ggplot(inset_long, aes(x = level, y = pct, fill = level)) +
  geom_col(width = 0.55, colour = "grey92", linewidth = 0.2) +
  geom_text(
    data = inset_sum,
    aes(x = level, y = pct_generalist, label = sprintf("%.1f%%", pct_generalist)),
    inherit.aes = FALSE,
    vjust = -0.2,
    size = 2.8,
    colour = "black"
  ) +
  scale_fill_manual(values = level_colors, guide = "none") +
  scale_y_continuous(limits = c(0, 105), breaks = c(0, 50, 100), labels = function(z) paste0(z, "%")) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7),
    axis.ticks = element_line(linewidth = 0.25),
    axis.line = element_line(linewidth = 0.25),
    plot.margin = margin(0, 0, 0, 0),
    panel.grid = element_blank()
  )

g_inset <- ggplotGrob(p_inset)
p_line_breadth <- p_line_breadth +
  annotation_custom(g_inset, xmin = 12.2, xmax = 19.7, ymin = 26, ymax = 75)

# Host phage load distribution: three lines, x = vOTU count (1..>=20), y = pct (%)
rinfection <- rbindlist(lapply(names(host_infection_dist_tables), function(lev) {
  dt <- copy(host_infection_dist_tables[[lev]])
  dt[, level := factor(paste0(toupper(substring(lev, 1, 1)), substring(lev, 2)), levels = c("Phylum", "Genus", "Species"))]
  dt[, x := as.integer(n_vOTUs)]
  dt[, y := n_hosts]
  dt[, .(level, x, y)]
}))
rinfection[, y_pct := 100 * y / sum(y), by = level]
p_line_infection <- ggplot(rinfection, aes(x = x, y = y_pct, colour = level)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_colour_manual(values = level_colors, name = "Level") +
  scale_x_continuous(breaks = x_breaks, labels = x_breaks_lab) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), labels = function(z) paste0(z, "%")) +
  labs(x = "Number of vOTUs per host", y = "Proportion (%)") +
  line_theme_pub()

ggsave("host_range_breadth_line_phylum_genus_species.pdf", p_line_breadth, width = 6, height = 4.5)
ggsave("host_range_breadth_line_phylum_genus_species.png", p_line_breadth, width = 6, height = 4.5, dpi = 150)
ggsave("host_infection_load_line_phylum_genus_species.pdf", p_line_infection, width = 6, height = 4.5)
ggsave("host_infection_load_line_phylum_genus_species.png", p_line_infection, width = 6, height = 4.5, dpi = 150)

# Combined figure (optional)
if (requireNamespace("gridExtra", quietly = TRUE)) {
  for (lev in c("phylum", "genus", "species")) {
    if (!is.null(host_range_plots[[lev]]) && !is.null(host_infection_plots[[lev]])) {
      p_comb <- gridExtra::grid.arrange(host_range_plots[[lev]], host_infection_plots[[lev]], ncol = 1)
      ggsave(paste0("host_range_and_infection_", lev, "_combined.pdf"), p_comb, width = 7, height = 9)
    }
  }
}
