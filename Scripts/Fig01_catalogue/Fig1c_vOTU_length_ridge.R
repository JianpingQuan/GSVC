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


input_file <- "file.path(PHAGE_ROOT, "finalFigures", "Figure1")/high_medium_quality_length.txt"
vOTU_length <- read.table(input_file, header = FALSE)
names(vOTU_length) <- c("Length", "Quality")

vOTU_length$Length <- as.numeric(vOTU_length$Length)
vOTU_length <- vOTU_length[complete.cases(vOTU_length), ]

# ======================
# 2. Preprocessing and colors
# ======================
library(ggplot2)
library(ggridges)
library(dplyr)

quality_levels <- c("Complete", "High-quality", "Medium-quality")
vOTU_length$Quality <- factor(vOTU_length$Quality, levels = quality_levels)

quality_colors <- c(
  "Medium-quality" = "#807DBA", 
  "High-quality" = "#4292C6",   
  "Complete" = "#41AB5D"         
)

# Compute 95% confidence interval
range_stats <- vOTU_length %>%
  group_by(Quality) %>%
  summarise(
    Lower = quantile(Length, 0.025),
    Upper = quantile(Length, 0.975),
    .groups = 'drop'
  ) %>%
  mutate(Quality_num = as.numeric(Quality))

# ======================
# 3. Plotting
# ======================
# ======================
# 3. Plotting (X-axis labels revised)
# ======================
# ======================
# 3. Plotting (log10 x-axis, kb text labels)
# ======================
p_combined <- ggplot(vOTU_length, aes(x = log10(Length), y = Quality, fill = Quality)) +
  
  # Ridge density layers
  stat_density_ridges(
    geom = "density_ridges_gradient",
    scale = 0.8, 
    alpha = 0.8,
    color = "black",              
    linewidth = 0.6,
    rel_min_height = 0.005,
    quantile_lines = TRUE,
    quantiles = c(0.25, 0.5, 0.75),
    calc_ecdf = TRUE,
    vline_color = "white",        
    vline_linetype = "dashed",    
    vline_width = 0.4             
  ) +
  
  # Red 95% range line
  geom_segment(
    data = range_stats,
    aes(x = log10(Lower), xend = log10(Upper), 
        y = Quality_num - 0.05, yend = Quality_num - 0.05),
    color = "#D73027", linewidth = 0.8, inherit.aes = FALSE 
  ) +
  geom_point(data = range_stats, aes(x = log10(Lower), y = Quality_num - 0.05),
             color = "#D73027", size = 1.2, inherit.aes = FALSE) +
  geom_point(data = range_stats, aes(x = log10(Upper), y = Quality_num - 0.05),
             color = "#D73027", size = 1.2, inherit.aes = FALSE) +
  
  # Key change: stats text kept as xx kb
  geom_text(
    data = range_stats,
    aes(x = log10(Lower) - 0.05, y = Quality_num - 0.08, 
        label = paste0(round(Lower/1000, 1), "-", round(Upper/1000, 1), " kb")),
    hjust = 1, vjust = 1, size = 3.5, fontface = "bold.italic", color = "grey20", inherit.aes = FALSE
  ) +
  
  scale_fill_manual(values = quality_colors) +
  
  # Key change: x-axis shows raw log10 values
  scale_x_continuous(
    name = expression(bold(log[10] * " Genome Length (bp)")),
    breaks = seq(3, 7, by = 1), 
    labels = seq(3, 7, by = 1),
    expand = expansion(mult = c(0.15, 0.05)) 
  ) +
  
  scale_y_discrete(expand = expansion(add = c(0.2, 0.2))) + 
  
  theme_bw() + 
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30"),
    plot.title.position = "plot", 
    panel.border = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    # X-axis grid lines
    panel.grid.major.x = element_line(color = "grey92", linetype = "dashed", linewidth = 0.3),
    axis.line.x = element_line(color = "black", linewidth = 0.6),
    axis.text.y = element_text(size = 15, face = "bold", color = "black"),
    axis.text.x = element_text(size = 12, color = "black"), 
    axis.ticks.y = element_blank(),
    legend.position = "none",
    plot.margin = margin(15, 15, 15, 15)
  ) +
  labs(
    title = "Distribution of vOTU Lengths",
    subtitle = "White dashed: Quartiles | Red bar: 95% interval in kb"
  )

# ======================
# 4. Save PDF
# ======================
output_path <- file.path(dirname(input_file), "high_medium_quality_length_distribution.pdf")

ggsave(
  filename = output_path, 
  plot = p_combined, 
  width = 9, # Slightly wider for more labels
  height = 5.5, 
  units = "in", 
  device = cairo_pdf
)
