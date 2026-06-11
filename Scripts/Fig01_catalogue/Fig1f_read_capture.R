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
library(tidyr)
library(dplyr)
library(ggsci)

# setwd: using PHAGE_ROOT above

# 1. Build raw data
raw_data <- data.frame(
  GSVC = c(15.93, 45.05, 29.28, 57.94, 45.75, 62.88, 52.32, 58.2, 64.15, 37.82, 56.04, 48.5, 38.7, 94.32, 10.88, 10.08),
  Study1 = c(5.02, 6.42, 4.38, 16.25, 14.58, 10.18, 13.02, 12.87, 20.27, 20.49, 14.26, 4.32, 13.23, 0.26, 2.4, 0.47),
  Study2 = c(11.16, 20.75, 8.92, 24.65, 23.41, 29.21, 20.89, 24.65, 27.13, 26.65, 17.64, 16.01, 14.7, 0.27, 2.82, 4.26),
  ICTV = c(0.52, 1.25, 1.26, 5.79, 2.7, 3.83, 3.39, 7.78, 5.53, 12.93, 2.57, 1.37, 1.22, 0.6, 0.55, 0.07),
  IMG_VR = c(28.91, 52.56, 37, 5.79, 60.1, 74.81, 69.08, 79.95, 80.21, 12.93, 78.99, 58.56, 47.22, 95.09, 13.48, 14.74)
)

# 2. Reshape data
df_long <- raw_data %>%
  pivot_longer(cols = everything(), names_to = "Database", values_to = "Percentage")

# 3. Set y-axis order
# ggplot draws y-axis bottom-up; top item goes last in levels
target_order <- rev(c("ICTV", "IMG_VR", "GSVC", "Study1", "Study2"))

# Factor Database with target order
df_long$Database <- factor(df_long$Database, levels = target_order)

# Order top to bottom: Study2, Study1, GSVC, IMG_VR, ICTV
my_colors <- c("#EEC78A", "#EEE9A2", "#CBE4B1", "#B3DDCB", "#B8E5FA")


# 4. Mean and standard deviation
df_summary <- df_long %>%
  group_by(Database) %>%
  summarise(
    mean_val = mean(Percentage),
    sd_val = sd(Percentage)
  )

# 5. Plot
p <- ggplot() +
  # Horizontal bars
  geom_bar(data = df_summary, aes(y = Database, x = mean_val, fill = Database), 
           stat = "identity", width = 0.7, alpha = 0.7, color = "black", linewidth = 0.3) +
  
  # Raw data points (gray hollow circles)
  geom_jitter(data = df_long, aes(y = Database, x = Percentage), 
              color = "grey40", shape = 1, size = 2, height = 0.15, alpha = 0.5) +
  
  # Horizontal error bars
  geom_errorbarh(data = df_summary, aes(y = Database, xmin = mean_val - sd_val, xmax = mean_val + sd_val), 
                 height = 0.2, color = "black", linewidth = 0.7) +
  
  # JCO-style palette
  scale_fill_manual(values = my_colors) +
  
  # Labels
  labs(
    x = "Alignment Percentage (%)",
    y = "Database Source",
    title = "Comparative Analysis of Reads Mapping"
  ) +
  
  # Theme
  theme_bw() + 
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 11, color = "black"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

# Display plot
print(p)

# PDF is vector format; preferred for submission
ggsave(
  filename = "Database_Comparison_Clean.pdf", 
  plot = p, 
  width = 7,         # Width
  height = 5,        # Height
  units = "in",      # Units: inches
  device = "pdf",    # Device type
  bg = "white"       # White background (avoids transparency in some viewers)
)
