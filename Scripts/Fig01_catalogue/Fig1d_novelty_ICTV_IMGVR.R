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




# 1. Load required packages
library(ggplot2)
library(tidyr)
library(dplyr)

# setwd: using PHAGE_ROOT above

# 2. Prepare data
# Values per description:
# Nucleotide: ICTV (0.97%), IMG/VR (21.6%)
# Protein: ICTV (11.5%), IMG/VR (46.9%)
data <- data.frame(
  Level = rep(c("Nucleotide Level", "Protein Level"), each = 2),
  Database = rep(c("ICTV", "IMG/VR"), 2),
  Percentage = c(0.97, 21.6, 11.5, 23.6)
)

# Lock facet order: Nucleotide left, Protein right
data$Level <- factor(data$Level, levels = c("Nucleotide Level", "Protein Level"))

# 3. Plot
p <- ggplot(data, aes(x = Database, y = Percentage, fill = Database)) +
  # Bar chart
  geom_bar(stat = "identity", width = 0.6, color = "black", size = 0.3) +
  # Facets: nucleotide left, protein right
  facet_wrap(~Level, scales = "fixed") + 
  # Percentage labels above bars
  geom_text(aes(label = paste0(Percentage, "%")), vjust = -0.5, size = 6) +
  # Academic blue and purple palette
  scale_fill_manual(values = c("ICTV" = "#7DAEE0", "IMG/VR" = "#B395BD")) +
  # Axis labels
  labs(x = NULL, y = "Percentage of GSVC Sequences (%)") +
  # Y-axis range with headroom for annotations
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  # Clean theme
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey90", color = NA), # Facet title background
    strip.text = element_text(size = 14),           # Facet title text
    axis.text = element_text(color = "black", size = 10),
    legend.position = "none",                                     # Hide legend (x-axis has labels)
    panel.spacing = unit(2, "lines"),  # Wider gap between facets
    axis.title.y = element_text(size = 16, color = "black", margin = margin(r = 10)), # Y-axis title
    axis.text.y = element_text(size = 14, color = "black"), # Y-axis ticks
    axis.text.x = element_text(size = 14, color = "black"), # X-axis labels (ICTV, etc.)
  )

# 4. Add large "Novel vOTUs" annotation
# Place above the Nucleotide Level facet
ann_text <- data.frame(
  Level = factor(c("Nucleotide Level", "Protein Level"), 
                 levels = c("Nucleotide Level", "Protein Level")),
  label = c("78.4% Novel vOTUs", "74.9% Novel Genes"), # Derived from data
  x = 1.5, 
  y = 92
)

# 2. Revised plot with annotations
p_final <- p + 
  # Horizontal line in each facet
  geom_segment(data = ann_text, 
               aes(x = 0.6, xend = 2.4, y = 85, yend = 85), 
               inherit.aes = FALSE, colour = "#299D8F", size = 0.5) +
  # Annotation text in each facet
  geom_text(data = ann_text, 
            aes(x = x, y = y, label = label), 
            inherit.aes = FALSE, color = "black", size = 5)

# Print result
print(p_final)

# 5. Save figure
ggsave("GSVC_Novelty_Plot.pdf", p_final, width = 6, height = 6)
