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


## 5) Count samples per country
# setwd: using PHAGE_ROOT above


# ========== 0. Load R packages ==========
# Load required packages
library(tidyverse)
library(maps)
library(ggrepel)
library(sf)
library(viridis)


# ========== Read China province boundaries ========== 
china_province <- st_read("China_sheng.geojson", quiet = TRUE)


# ========== Sample count data ==========
df <- data.frame(
  Country_clean = c("Australia","Austria","Belgium","Brazil", "Bulgaria","Canada", "China", "Denmark",
                    "France", "Germany","Ghana", "Ireland","Italy","Netherlands","New Zealand", 
                    "Norway","Poland", "South Korea", "Spain", "Thailand", "USA", 
                    "United Kingdom", "Hungary", "Russia", "Japan"),
  Sample_Count = c(830,39,22,86,21,305, 1811, 449, 432, 44,30, 89,20,21,45, 29,20, 90, 361, 116, 308, 82, 13, 16, 60)
)


# World map (excluding Antarctica)
world_map <- map_data("world") %>% filter(region != "Antarctica")


# Country name mapping (assign Taiwan to China)
df <- df %>%
  mutate(Country_map = case_when(
    Country_clean == "USA" ~ "USA",
    Country_clean == "United Kingdom" ~ "UK",
    Country_clean == "South Korea" ~ "South Korea",
    Country_clean == "New Zealand" ~ "New Zealand",
    Country_clean == "Ghana" ~ "Ghana",
    Country_clean == "China" ~ "China",  # Explicit China label
    TRUE ~ Country_clean
  ))


# Map Taiwan region to China in world map data
world_map <- map_data("world") %>%
  mutate(region = ifelse(region == "Taiwan", "China", region)) %>%
  filter(region != "Antarctica")


# Merge map with sample counts
map_data_colored <- world_map %>% left_join(df, by = c("region" = "Country_map"))

# ========== Read sampling locations ========== 
loc_raw <- readLines("locations.txt")
loc_df <- data.frame(raw = loc_raw) %>%
  separate(raw, into = c("Country", "City", "Coord"), sep = ":", fill = "right", extra = "merge") %>%
  separate(Coord, into = c("lat", "lat_dir", "long", "long_dir"), sep = " ") %>%
  mutate(
    lat = as.numeric(lat) * ifelse(lat_dir == "S", -1, 1),
    long = as.numeric(long) * ifelse(long_dir == "W", -1, 1)
  ) %>%
  select(Country, City, lat, long)

# ========== Unified coloring (including China) ==========
# Replace region with China in map data; sample counts already mapped
china_map <- china_province %>%
  mutate(region = "China") %>%
  left_join(df, by = c("region" = "Country_clean"))

# ========== Cool-tone color scheme ==========
# Scheme: Deep Sea blues
# From pale teal (no/low data) to teal, then deep blue (high values)
cool_colors <- c("#e0f2f1", "#80deea", "#26c6da", "#00838f", "#00363a")

# ========== Plotting ========== 
p <- ggplot() +
  # Ocean background (optional polish)
  geom_rect(aes(xmin = -180, xmax = 180, ymin = -60, ymax = 85), 
            fill = "#f8f9fa", color = NA) + 
  
  # World base map
  geom_polygon(
    data = map_data_colored,
    aes(x = long, y = lat, group = group, fill = Sample_Count),
    color = "white", linewidth = 0.1 # White borders for a cleaner look
  ) +
  
  # Overlay China province boundaries
  geom_sf(data = china_province, fill = NA, color = "#455a64", linewidth = 0.2) +
  
  # Sampling points: high-contrast accent color
  geom_point(
    data = loc_df,
    aes(x = long, y = lat),
    color = "#ff8f00", size = 0.6, alpha = 0.7 # Amber accent on cool base map
  ) +
  
  # Cool-tone color scale
  scale_fill_gradientn(
    colours = cool_colors,
    na.value = "#eceff1", # Light gray-blue for no-data regions
    name = "Sample Count",
    trans = "log10", # log10 scale for more even color distribution
    guide = guide_colorbar(
      barwidth = 1, barheight = 10,
      frame.colour = "black", ticks.colour = "black"
    )
  ) +
  
  # Projection and axes
  coord_sf(
    crs = "EPSG:4326",
    xlim = c(-180, 180),
    ylim = c(-60, 85),
    expand = FALSE
  ) +
  
  # Theme: minimal style
  theme_minimal(base_family = "sans") +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 9, color = "#546e7a"),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "#f8f9fa", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    
    # Legend position
    legend.position = c(0.12, 0.3), # Place legend in lower-left whitespace
    legend.background = element_rect(fill = alpha("white", 0.7)),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 8),
    
    # Title aligned left
    plot.title = element_text(face = "bold", size = 14, color = "#263238", hjust = 0),
    plot.margin = margin(10, 10, 10, 10)
  )

# Custom tick marks (colors fine-tuned)
p <- p + 
  geom_segment(data = data.frame(lat = lat_ticks),
               aes(x = -180, xend = -176, y = lat, yend = lat),
               color = "#90a4ae") +
  geom_segment(data = data.frame(long = long_ticks),
               aes(x = long, xend = long, y = 85, yend = 81),
               color = "#90a4ae")

print(p)

# Save figure
ggsave("country_sample_counts_by_continent_positions.pdf", p, width = 8, height = 6, units = "in", dpi = 600)