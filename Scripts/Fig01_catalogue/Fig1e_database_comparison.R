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




library(ggplot2)
library(showtext)
# setwd: using PHAGE_ROOT above

# 1. Configure font (Arial)
font_add_google("Arimo", "arial")
showtext_auto()

# 2. Prepare data
# Reshape to long format for ggplot2
data <- data.frame(
  Study = rep(c("GSVC", "Study1", "Study2"), each = 3),
  Quality = rep(c("Medium-quality", "High-quality", "Complete"), 3),
  Count = c(115842, 54486, 54277,   # GSVC
            19623, 7747, 5149,      # Study1
            26184, 9600, 12515)     # Study2
)

# 3. Lock factor order
# X-axis study order
data$Study <- factor(data$Study, levels = c("GSVC", "Study1", "Study2"))
# Stack order bottom to top: Medium -> High -> Complete
data$Quality <- factor(data$Quality, levels = c("Medium-quality", "High-quality", "Complete"))

# 4. Plot
p <- ggplot(data, aes(x = Study, y = Count, fill = Quality)) +
  # Stacked bars with thin black borders
  geom_bar(stat = "identity", width = 0.7, color = "black", size = 0.3) +
  # Publication-quality color palette
  scale_fill_manual(values = c("Medium-quality" = "#A8DADC", 
                               "High-quality" = "#457B9D", 
                               "Complete" = "#1D3557")) +
  # Axes and title
  labs(x = NULL, y = "Number of vOTUs", fill = "Genome Quality") +
  # Clean theme
  theme_classic() +
  theme(
    text = element_text(family = "arial"),
    # Larger axis text
    axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 14, color = "black"),
    # Legend styling
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "top"
  )

# 5. Save figure
ggsave("vOTU_Quality_Stacked_Bar.pdf", p, width = 7, height = 8, device = cairo_pdf)

# Display
print(p)
