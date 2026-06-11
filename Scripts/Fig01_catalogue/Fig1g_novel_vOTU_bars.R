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



# Load libraries
library(ggplot2)
library(dplyr)

# 1. Prepare data
data <- data.frame(
  category = c("Database-matched", "PreviousStudy-matched", "GSVC-novel"),
  count = c(48575, 29972, 148791) 
)

# 2. Compute percentages and label positions (key step)
data <- data %>%
  mutate(category = factor(category, levels = c("Database-matched", "PreviousStudy-matched", "GSVC-novel"))) %>%
  arrange(desc(category)) %>% # Sort to match stack order
  mutate(
    percentage = count / sum(count),
    # Midpoint of each slice for label placement
    pos = cumsum(count) - (count / 2),
    label_text = paste0(format(count, big.mark=","), "\n(", round(percentage * 100, 1), "%)")
  )

# 3. Set colors
my_colors <- c("Database-matched" = "#74AED4", 
               "PreviousStudy-matched" = "#CFAFD4", 
               "GSVC-novel" = "#D3E2B7") 

# 4. Plot
p <- ggplot(data, aes(x = 3, y = count, fill = category)) +
  geom_bar(stat = "identity", color = "white", width = 1) +
  coord_polar(theta = "y") +
  # Text label layer
  geom_text(aes(y = pos, label = label_text), 
            x = 3,           # x=3 centers text in each ring segment
            color = "black", # Adjust for background contrast if needed
            size = 4.5) +
  xlim(0.8, 3.5) +
  theme_void() + 
  scale_fill_manual(values = my_colors) +
  annotate("text", x = 0.8, y = 0, label = paste0("Total vOTUs\n", "224,605"), 
           size = 6, fontface = "bold") +
  labs(title = "Novelty breakdown of the GSVC",
       fill = "Category") +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.title = element_text(face = "bold"))

# 5. Save
ggsave("Figure_1g_Novelty_Breakdown_with_labels.pdf", 
       plot = p, 
       width = 7,      # Extra width for legend
       height = 6, 
       device = cairo_pdf) # cairo_pdf embeds fonts better
