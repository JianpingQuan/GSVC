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


#!/usr/bin/env Rscript
# Multi-family polar phylogeny — R version (ape + ggplot2), same approach as plot_multi_family_circular_phylo.py
# Dependencies: ape, readr, dplyr, ggplot2, stringr, grid
#
# Rscript scripts/plot_multi_family_circular_phylo.R --root path/to/vcontact3_run2
# Rscript ... --families "Microviridae,Inoviridae" --out path/prefix_no_ext --full_circle
# Rscript ... --no_sector_fill   # Disable pale sector fill background

suppressPackageStartupMessages({
  library(ape)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(grid)
})

COLOR_NOVEL <- "#2CA25A"
COLOR_KNOWN <- "#000000"
COLOR_INTERNAL <- "#000000"
R_TREE_OUTER <- 0.92

# Must match actual *family*_tree.nw files under exports/newick; this run is mostly Caudoviricetes etc.; when Adeno/Calici… are absent, commonly use e.g.:
DEFAULT_FAMILIES <- c(
  "Microviridae", "Inoviridae", "Zobellviridae",
  "Tectiviridae", "Herelleviridae", "Autographiviridae"
)

strip_nhx <- function(txt) {
  str_replace_all(txt, "\\[&&NHX:[^\\]]*\\]", "")
}

is_ref_tip <- function(name) {
  n <- trimws(name)
  startsWith(n, "NC_") || startsWith(n, "NZ_") || startsWith(n, "NR_") ||
    startsWith(n, "XM_") || startsWith(n, "NW_") || startsWith(n, "NG_")
}

norm_ref <- function(s) {
  tolower(trimws(as.character(s))) %in% c("true", "1", "yes")
}

novel_or_unplaced <- function(val) {
  s <- tolower(as.character(val))
  grepl("novel", s) || grepl("unplaced", s)
}

load_assignments <- function(path) {
  d <- read_csv(path, show_col_types = FALSE)
  d <- d[!duplicated(d$Genome), , drop = FALSE]
  d <- as.data.frame(d)
  rownames(d) <- as.character(d$Genome)
  d
}

is_novel_tip <- function(tip_name, adf) {
  if (is_ref_tip(tip_name)) {
    return(FALSE)
  }
  if (!tip_name %in% rownames(adf)) {
    return(TRUE)
  }
  row <- adf[tip_name, , drop = FALSE]
  if ("Reference" %in% names(row) && norm_ref(row$Reference[1])) {
    return(FALSE)
  }
  for (cn in c(
    "class_prediction", "order_prediction", "family_prediction",
    "subfamily_prediction", "genus_prediction"
  )) {
    if (cn %in% names(row) && novel_or_unplaced(row[[cn]][1])) {
      return(TRUE)
    }
  }
  FALSE
}

find_family_newick <- function(nd, family) {
  exact <- file.path(nd, paste0(family, "_tree.nw"))
  if (file.exists(exact)) {
    return(exact)
  }
  fl <- list.files(nd, pattern = "tree\\.nw$", full.names = TRUE)
  fam_l <- tolower(family)
  hit <- fl[grepl(fam_l, tolower(basename(fl)), fixed = TRUE)]
  if (length(hit)) {
    return(hit[which.min(nchar(basename(hit)))[1]])
  }
  character(0)
}

children_ordered <- function(tr) {
  mx <- max(tr$edge)
  ch <- vector("list", mx)
  for (i in seq_len(nrow(tr$edge))) {
    p <- tr$edge[i, 1L]
    c <- tr$edge[i, 2L]
    ch[[p]] <- c(ch[[p]], c)
  }
  ch
}

root_node <- function(tr) {
  setdiff(unique(tr$edge[, 1L]), tr$edge[, 2L])[1L]
}

to_nested <- function(tr, node, chl) {
  n <- Ntip(tr)
  if (node <= n) {
    return(list(tip = TRUE, node = node, name = tr$tip.label[node]))
  }
  list(
    tip = FALSE, node = node,
    children = lapply(chl[[node]], function(k) to_nested(tr, k, chl))
  )
}

terminal_labels_inorder <- function(clade) {
  if (clade$tip) {
    return(clade$name)
  }
  unlist(lapply(clade$children, terminal_labels_inorder), use.names = FALSE)
}

sector_palette <- function(n) {
  if (n <= 0L) {
    return(character(0))
  }
  hcl.colors(max(n, 2L), palette = "Dark 3")[seq_len(n)]
}

read_tree_clean <- function(path) {
  raw <- paste(readLines(path, warn = FALSE), collapse = "")
  raw <- strip_nhx(raw)
  tr <- read.tree(text = raw)
  ladderize(tr, right = TRUE)
}

#' Compute node y (as in Python: assign order to tips first; internal nodes get mean y of first and last child)
compute_y_node <- function(tr, nest) {
  n <- Ntip(tr)
  rt <- root_node(tr)
  chl <- children_ordered(tr)
  mx <- max(tr$edge)
  yv <- numeric(mx)
  tips_rev <- rev(terminal_labels_inorder(nest))
  for (i in seq_along(tips_rev)) {
    ti <- which(tr$tip.label == tips_rev[i])[1L]
    yv[ti] <- length(tips_rev) - i + 1L
  }
  post_int <- function(node) {
    if (node <= n) {
      return(yv[node])
    }
    kids <- chl[[node]]
    if (is.null(kids) || !length(kids)) {
      return(NA_real_)
    }
    for (k in kids) {
      post_int(k)
    }
    yv[node] <<- (yv[kids[1L]] + yv[kids[length(kids)]]) / 2
    yv[node]
  }
  post_int(rt)
  yv
}

collect_wedge_geom <- function(tr, adf, theta_lo, theta_hi, sector_hex, lw) {
  n <- Ntip(tr)
  chl <- children_ordered(tr)
  rt <- root_node(tr)
  nest <- to_nested(tr, rt, chl)
  y_node <- compute_y_node(tr, nest)

  nnode <- Ntip(tr) + tr$Nnode
  xd <- node.depth.edgelength(tr)
  if (is.null(xd) || !length(xd) || max(xd, na.rm = TRUE) <= 0) {
    xd <- node.depth(tr)
  }
  if (length(xd) < nnode) {
    xd <- c(as.numeric(xd), rep(0, nnode - length(xd)))
  }
  xd <- as.numeric(xd)

  tip_idx <- seq_len(n)
  ymin <- min(y_node[tip_idx]) - 0.5
  ymax <- max(y_node[tip_idx]) + 0.5
  yspan <- ymax - ymin
  if (yspan <= 0) {
    yspan <- 1
  }

  xmax <- max(xd, na.rm = TRUE)
  if (xmax <= 0) {
    xmax <- 1
  }

  margin_ang <- 0.035 * (theta_hi - theta_lo)
  th_lo <- theta_lo + margin_ang
  th_hi <- theta_hi - margin_ang

  theta_of_y <- function(y) {
    th_lo + (th_hi - th_lo) * (y - ymin) / yspan
  }
  r_of_x <- function(x) {
    (x / xmax) * R_TREE_OUTER
  }

  segs <- list()
  arcs <- list()
  aid <- 1L

  draw_clade <- function(clade, x_start) {
    node <- clade$node
    x_here <- xd[node]
    y_here <- y_node[node]
    th <- theta_of_y(y_here)
    r_s <- r_of_x(x_start)
    r_e <- r_of_x(x_here)
    col <- if (clade$tip) {
      if (is_novel_tip(clade$name, adf)) COLOR_NOVEL else COLOR_KNOWN
    } else {
      COLOR_INTERNAL
    }
    segs[[length(segs) + 1L]] <<- tibble(th = th, r0 = r_s, r1 = r_e, col = col, lw = lw)

    if (!clade$tip && length(clade$children)) {
      y_top <- y_node[clade$children[[1L]]$node]
      y_bot <- y_node[clade$children[[length(clade$children)]]$node]
      t_a <- theta_of_y(y_bot)
      t_b <- theta_of_y(y_top)
      t_lo <- min(t_a, t_b)
      t_hi <- max(t_a, t_b)
      nn <- max(8L, as.integer(32 * abs(t_hi - t_lo) / pi))
      arc_th <- seq(t_lo, t_hi, length.out = nn)
      arcs[[length(arcs) + 1L]] <<- tibble(
        th = arc_th, r = r_of_x(x_here), col = COLOR_INTERNAL, lw = lw, gid = aid
      )
      aid <<- aid + 1L
      for (ch in clade$children) {
        draw_clade(ch, x_here)
      }
    }
  }

  draw_clade(nest, 0)

  seg_df <- if (length(segs)) bind_rows(segs) else tibble()
  arc_df <- if (length(arcs)) bind_rows(arcs) else tibble()

  ths <- vapply(tip_idx, function(ti) theta_of_y(y_node[ti]), numeric(1))
  labs <- tr$tip.label
  si <- order(ths)
  dth <- diff(ths[si])
  wid <- if (length(dth) == 0) {
    rep(0.02, length(ths))
  } else {
    pmax(pmin(c(dth[1], dth), c(dth, dth[length(dth)])) * 0.42, 0.0025)
  }
  w_by_orig <- numeric(length(ths))
  for (j in seq_along(si)) {
    w_by_orig[si[j]] <- wid[j]
  }
  ring_df <- tibble(
    th = ths,
    lab = labs,
    w = w_by_orig,
    fill = sector_hex
  ) %>%
    mutate(
      r0 = R_TREE_OUTER + 0.02,
      r1 = r0 + 0.07,
      xmin = th - w / 2,
      xmax = th + w / 2
    )

  tips_df <- tibble(
    th = ths,
    r = vapply(tip_idx, function(ti) r_of_x(xd[ti]), numeric(1)),
    fill = sector_hex
  )

  list(seg = seg_df, arc = arc_df, tips = tips_df, ring = ring_df, n_tips = n)
}

draw_depth_arcs_gg <- function(theta_lo, theta_hi, arc_lw) {
  th <- seq(theta_lo, theta_hi, length.out = max(64L, as.integer(120 * (theta_hi - theta_lo) / pi)))
  bind_rows(lapply(c(0.25, 0.5, 0.75, 1.0), function(frac) {
    tibble(th = th, r = R_TREE_OUTER * frac, lw = arc_lw, ring_id = frac)
  }))
}

# ggplot2::coord_polar treats geom_segment as planar segments then transforms them; radial branches become chords through center → messy branches.
# Same as matplotlib: draw radial segments and arcs in Cartesian plane with x=r*cos(theta), y=r*sin(theta).
polar_seg_to_cart <- function(d) {
  if (!nrow(d)) {
    return(d)
  }
  d %>%
    mutate(
      x = r0 * cos(th), y = r0 * sin(th),
      xend = r1 * cos(th), yend = r1 * sin(th)
    )
}

polar_path_to_cart <- function(d, thcol = "th", rcol = "r", groupcol = "gid") {
  if (!nrow(d)) {
    return(tibble(x = numeric(0), y = numeric(0), gid = numeric(0), col = character(0), lw = numeric(0)))
  }
  parts <- split(seq_len(nrow(d)), d[[groupcol]])
  bind_rows(lapply(parts, function(ii) {
    block <- d[ii, , drop = FALSE]
    tibble(
      x = block[[rcol]] * cos(block[[thcol]]),
      y = block[[rcol]] * sin(block[[thcol]]),
      gid = block[[groupcol]][1],
      col = block$col[1],
      lw = block$lw[1]
    )
  }))
}

polar_guide_to_cart <- function(d) {
  if (!nrow(d)) {
    return(tibble(x = numeric(0), y = numeric(0), guide_gid = character(0), lw = numeric(0)))
  }
  d %>%
    mutate(guide_gid = as.character(interaction(ring_id, r, drop = TRUE))) %>%
    group_by(guide_gid) %>%
    reframe(x = r * cos(th), y = r * sin(th), lw = lw[1])
}

ring_to_polygons <- function(ring_df) {
  if (!nrow(ring_df)) {
    return(tibble(x = numeric(0), y = numeric(0), pid = integer(0), fill = character(0)))
  }
  ring_df %>%
    mutate(pid = row_number()) %>%
    group_by(pid, fill) %>%
    reframe(
      x = c(
        r0[1] * cos(th[1] - w[1] / 2), r1[1] * cos(th[1] - w[1] / 2),
        r1[1] * cos(th[1] + w[1] / 2), r0[1] * cos(th[1] + w[1] / 2)
      ),
      y = c(
        r0[1] * sin(th[1] - w[1] / 2), r1[1] * sin(th[1] - w[1] / 2),
        r1[1] * sin(th[1] + w[1] / 2), r0[1] * sin(th[1] + w[1] / 2)
      )
    )
}

tips_to_cart <- function(d) {
  if (!nrow(d)) {
    return(d)
  }
  d %>% mutate(x = r * cos(th), y = r * sin(th))
}

# Pale sector background (reference semicircle multi-panel plots with per-group fill)
sector_wedge_df <- function(sid, th0, th1, r_outer, fill_base, alpha = 0.14, n = 56) {
  th_seq <- seq(th0, th1, length.out = n)
  fc <- grDevices::adjustcolor(fill_base, alpha.f = alpha)
  tibble(
    sid = sid,
    x = c(0, r_outer * cos(th_seq), 0),
    y = c(0, r_outer * sin(th_seq), 0),
    sfill = fc
  )
}

get_arg <- function(argv, name, default) {
  i <- match(name, argv)
  if (is.na(i) || i >= length(argv)) {
    return(default)
  }
  argv[[i + 1]]
}

argv <- commandArgs(trailingOnly = TRUE)
root <- get_arg(argv, "--root", NA_character_)
if (is.na(root)) {
  root <- ""/vcontact3_run2"
}
fams_raw <- get_arg(argv, "--families", paste(DEFAULT_FAMILIES, collapse = ","))
out_prefix <- get_arg(argv, "--out", NA_character_)
full_circle <- ("--full_circle" %in% argv) || ("--full-circle" %in% argv)
no_sector_fill <- ("--no_sector_fill" %in% argv) || ("--no-sector-fill" %in% argv)

nd <- file.path(root, "exports", "newick")
fa <- file.path(root, "exports", "final_assignments.csv")
if (is.na(out_prefix)) {
  out_prefix <- file.path(root, "postanalysis", "figures_nm", "Fig_multi_family_circular_phylo_R")
}

if (!file.exists(fa)) {
  stop("No final_assignments.csv: ", fa, call. = FALSE)
}

families <- strsplit(fams_raw, ",", fixed = TRUE)[[1]]
families <- trimws(families[nzchar(families)])

adf <- load_assignments(fa)
resolved <- list()
missing <- character(0)

for (fam in families) {
  p <- find_family_newick(nd, fam)
  if (length(p) == 0) {
    missing <- c(missing, fam)
    next
  }
  tr <- tryCatch(read_tree_clean(p), error = function(e) NULL)
  if (is.null(tr)) {
    missing <- c(missing, paste0(fam, " (read fail)"))
    next
  }
  nt <- Ntip(tr)
  if (nt == 0L) {
    missing <- c(missing, paste0(fam, " (empty)"))
    next
  }
  resolved[[length(resolved) + 1L]] <- list(family = fam, path = p, tree = tr, n = nt)
}

if (!length(resolved)) {
  stop("No match family newick", call. = FALSE)
}

n_total <- sum(vapply(resolved, function(x) x$n, 1L))
n_w <- length(resolved)
cols_sec <- sector_palette(n_w)

angle_span <- if (full_circle) 2 * pi else pi
gap <- min(0.04, 0.5 / max(n_total, 1)) * (angle_span / max(n_w, 1))
usable <- angle_span - gap * (n_w + 1)
theta <- gap

lw <- 0.35
arc_lw <- max(lw * 0.55, 0.28)

all_seg <- list()
all_arc <- list()
all_tips <- list()
all_ring <- list()
th0_list <- numeric(n_w)
th1_list <- numeric(n_w)

for (i in seq_along(resolved)) {
  item <- resolved[[i]]
  span <- usable * (item$n / n_total)
  th0 <- theta
  th1 <- theta + span
  theta <- th1 + gap
  th0_list[i] <- th0
  th1_list[i] <- th1

  g <- collect_wedge_geom(item$tree, adf, th0, th1, cols_sec[i], lw)
  all_seg[[i]] <- g$seg
  all_arc[[i]] <- g$arc
  all_tips[[i]] <- g$tips
  all_ring[[i]] <- g$ring
}

guide_df <- if (full_circle) {
  draw_depth_arcs_gg(0, angle_span, arc_lw)
} else {
  draw_depth_arcs_gg(gap, angle_span - gap, arc_lw)
}

seg_df <- bind_rows(all_seg)
arc_df <- bind_rows(all_arc)
tips_df <- bind_rows(all_tips)
ring_df <- bind_rows(all_ring)

r_max <- 1.02 + 0.07 + 0.1
lim <- r_max * 1.08

seg_c <- polar_seg_to_cart(seg_df)
arc_c <- polar_path_to_cart(arc_df, "th", "r", "gid")
guide_c <- polar_guide_to_cart(guide_df)
ring_poly <- ring_to_polygons(ring_df)
tips_c <- tips_to_cart(tips_df)

wedge_layers <- if (!no_sector_fill) {
  bind_rows(lapply(seq_len(n_w), function(i) {
    sector_wedge_df(i, th0_list[i], th1_list[i], lim * 0.998, cols_sec[i], alpha = 0.14)
  }))
} else {
  tibble()
}

p <- ggplot()

if (nrow(wedge_layers)) {
  p <- p + geom_polygon(
    data = wedge_layers,
    aes(x = x, y = y, group = sid, fill = sfill),
    colour = NA
  )
}

if (nrow(guide_c)) {
  p <- p + geom_path(
    data = guide_c,
    aes(x = x, y = y, group = guide_gid),
    colour = "#9e9e9e", linewidth = arc_lw * 0.35,
    linetype = "dashed", alpha = 0.85, lineend = "round"
  )
}

if (nrow(seg_c)) {
  p <- p + geom_segment(
    data = seg_c,
    aes(x = x, y = y, xend = xend, yend = yend, colour = col, linewidth = lw),
    lineend = "round"
  )
}

if (nrow(arc_c)) {
  p <- p + geom_path(
    data = arc_c,
    aes(x = x, y = y, group = gid, colour = col, linewidth = lw),
    lineend = "round"
  )
}

if (nrow(ring_poly)) {
  p <- p + geom_polygon(
    data = ring_poly,
    aes(x = x, y = y, group = pid, fill = fill),
    colour = NA
  )
}

if (nrow(tips_c)) {
  p <- p + geom_point(
    data = tips_c,
    aes(x = x, y = y, fill = fill),
    shape = 21, colour = "#333333", stroke = 0.18, size = 1.6
  )
}

p <- p +
  scale_colour_identity() +
  scale_fill_identity() +
  scale_linewidth_identity() +
  coord_fixed(
    ratio = 1,
    xlim = c(-lim, lim),
    ylim = if (full_circle) c(-lim, lim) else c(-0.05 * lim, lim),
    expand = FALSE,
    clip = "off"
  ) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    title = paste0(
      "Multi-family phylogeny (R Cartesian; ",
      if (full_circle) "full circle" else "semicircle",
      ") — green = novel/unplaced terminal"
    )
  ) +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5, margin = margin(b = 4)))

dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)

cap <- c(
  "Fig_multi_family_circular_phylo (R). vContact3 newick + final_assignments.",
  "Rendering: Cartesian x=r*cos(theta), y=r*sin(theta) (NOT ggplot2 coord_polar) to avoid spurious chord segments.",
  paste("Root:", root),
  paste("Families drawn:", paste(vapply(resolved, function(x) x$family, ""), collapse = ", ")),
  if (length(missing)) paste("Skipped:", paste(missing, collapse = "; ")) else NULL
)
writeLines(cap[!sapply(cap, is.null)], paste0(out_prefix, ".caption.txt"))

pdf_dev <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
ww <- if (full_circle) 11 else 12
hh <- if (full_circle) 11 else 6.6
ggsave(paste0(out_prefix, ".pdf"), p, width = ww, height = hh, device = pdf_dev, dpi = 300)
ggsave(paste0(out_prefix, ".png"), p, width = ww, height = hh, dpi = 200, bg = "white")

message("[OK] ", normalizePath(paste0(out_prefix, ".pdf"), winslash = "/"))
