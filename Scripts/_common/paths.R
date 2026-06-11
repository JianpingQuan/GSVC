# Shared paths for GSVC GitHub repository scripts.
args <- commandArgs(trailingOnly = FALSE)
sp <- sub("^--file=", "", args[grep("^--file=", args)])
if (length(sp) && nzchar(sp)) {
  n_up <- length(strsplit(dirname(sp), .Platform$file.sep)[[1]]) -
    length(strsplit(normalizePath(file.path(dirname(sp), "../../.."), winslash = "/"), "/")[[1]])
  # Default: script lives under Scripts/<subdir>/
  GSVC_ROOT <- normalizePath(file.path(dirname(sp), "..", ".."), winslash = "/", mustWork = FALSE)
} else {
  GSVC_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
for (ev in c("GSVC_ROOT", "PHAGE_PROJECT_ROOT")) {
  alt <- Sys.getenv(ev, "")
  if (nzchar(alt) && dir.exists(alt)) {
    GSVC_ROOT <- normalizePath(alt, winslash = "/")
    break
  }
}
if (!file.exists(file.path(GSVC_ROOT, "Pre-processed_Files/metadata/meta_augmented_combined.tsv"))) {
  warning("GSVC_ROOT may be wrong; set GSVC_ROOT or PHAGE_PROJECT_ROOT")
}
PHAGE_ROOT <- GSVC_ROOT
OUT_DIR <- file.path(GSVC_ROOT, "Scripts", "output")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(GSVC_ROOT)
