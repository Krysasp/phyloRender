############################################################
## SCRIPT CREDENTIALS
############################################################
## Tool       :  # phyloRender as part of the 'phyloAct' program
## Current script    :  # main.R
## Description       :  # This script tests the important features of MCC-tree plotting functionality
## Author            :  # Jonathan Chan
## Affiliation       :  # Institute of Health and Community Medicine
## Contact Email     :  # jonatasp92@gmail.com
## Date Created      :  # 2023-09-16
## Last Updated      :  # 2025-12-29
## Version           :  # v1.1-beta
## R Version         :  # Recommended on R version 4.4.2
## Required Packages :  # ape, ggtree, treeio, maps, phangorn, ggplot2
## Input Files       :  # requires mcc tr in nexus format
## Output Files      :  # plots, tables, animations
############################################################
############################################################
## main.R — Comprehensive MCC Tree Analysis & Animation
## Features: branch HPDs, rate/emergence labeling,
##           spatial plotting, discrete traits, animation,
##           BEAST1 & BEAST2 support, package/source-ready
############################################################

############################################################
## 0. SETUP ENVIRONMENT
############################################################

# Required CRAN packages
required_pkgs <- c("ape", "maps", "treeio", "phangorn", "ggplot2")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop("Missing required package: ", p)
  }
}
library(ape)
library(maps)
library(treeio)
library(phangorn)
library(ggplot2)

# Create output directories
FRAME_DIR <- "output/frames"
PLOT_DIR  <- "output/plots"
dir.create(FRAME_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(PLOT_DIR,  recursive=TRUE, showWarnings=FALSE)

############################################################
## 1. SOURCE FUNCTION SCRIPTS
############################################################

# Source all R scripts containing functions:
# - pts_and_hpds_at_time.R
# - plot_discrete_tree.R
# - plot_mcc_tree_with_hpds.R
# - latlon_hpds.R
# - discrete_traits.R
# These scripts include the following features:
#   1. Branch HPD extraction (emergence and rate)
#   2. Spatial HPDs & ancestral coordinates
#   3. Discrete trait extraction & colouring
#   4. Animation driver
#   5. BEAST1/BEAST2 hybrid reader
rfiles <- list.files("R", pattern="\\.R$", full.names=TRUE)
for (f in rfiles) {
  message("Sourcing: ", f)
  source(f)
}

############################################################
## 2. INPUT FILES & PARAMETERS
############################################################

# ===== USER EDITABLE PARAMETERS =====
# Path to BEAST MCC tree (BEAST1 or BEAST2)
MCC_FILE <- "data/example_mcc.tree"

# Select which discrete trait to visualize (integer index)
DISCRETE_TRAIT_INDEX <- 1

# Animation parameters
ANIMATION_FRAMES <- 30           # total number of frames
ANIMATION_PROP <- 1               # discrete trait prop index
ANIMATION_PREFIX <- "frame"       # prefix for saved PNG frames

# Branch HPD labeling options
BRANCH_LABEL_MODE <- c("emergence","rate") # modes to plot
NODES_TO_LABEL <- NULL   # if NULL, automatically pick internal nodes

# Spatial map limits (longitude, latitude)
MAP_XLIM <- c(-180, 180)
MAP_YLIM <- c(-60, 80)
# ====================================

############################################################
## 3. READ AND PREPARE MCC TREE
############################################################

# Unified MCC tree reader
# Features:
# - Reads BEAST1 or BEAST2 tree
# - Extracts lat/lon, HPDs, discrete traits, branch rates
tr <- read_beast_mcc(
  file     = MCC_FILE,
  spatial  = TRUE,   # extract lat/lon + HPDs
  discrete = TRUE,   # extract discrete traits
  rates    = TRUE    # extract substitution rates
)

# Optional: flip coordinates if necessary
# tr <- switch_latlon(tr)

# Calculate node times if missing
if(is.null(tr$nodeTimes)){
  tr$nodeTimes <- nodeHeights(tr) # decimal years
}

############################################################
## 4. PLOT BASIC TREE TO CHECK TOPOLOGY
############################################################

# Plot without tip labels
png(file.path(PLOT_DIR, "tree_topology.png"), 1200, 800)
plot.phylo(tr, show.tip.label=FALSE)
title("MCC Tree Topology")
dev.off()
# Use this to verify tree structure before plotting HPDs or traits

############################################################
## 5. PLOT DISCRETE TRAITS
############################################################

# Features:
# - Color tips and branches by discrete trait
# - Optionally add legend
# - `propIndex` selects which trait to display

png(file.path(PLOT_DIR, "tree_discrete_trait.png"), 1200, 800)
plot_discrete_tree(
  tr,
  propIndex = DISCRETE_TRAIT_INDEX,
  legpos    = "bottomleft",
  edgeCols  = TRUE
)
dev.off()

############################################################
## 6. PLOT SPATIAL MCC TREE WITH HPDs
############################################################

# Features:
# - Ancestral location points
# - HPD polygons around nodes
# - Optional discrete coloring
# - Map longitude/latitude limits adjustable

png(file.path(PLOT_DIR, "spatial_hpds.png"), 1200, 800)
plot_mcc_tree_with_hpds(
  tr,
  xlim = MAP_XLIM,
  ylim = MAP_YLIM,
  propIndex = DISCRETE_TRAIT_INDEX
)
title("Spatial MCC Tree with HPDs")
dev.off()

############################################################
## 7. BRANCH HPD LABELING
############################################################

# If NODES_TO_LABEL not specified, pick first 5 internal nodes
if(is.null(NODES_TO_LABEL)){
  NODES_TO_LABEL <- (length(tr$tip.label)+1):(length(tr$tip.label)+5)
}

# Loop over modes: emergence and substitution rate
for(mode in BRANCH_LABEL_MODE){
  png(file.path(PLOT_DIR, paste0("branch_", mode, "_hpds.png")), 1200, 800)
  plot_mcc_tree_with_branch_hpds(
    tr,
    nodes = NODES_TO_LABEL,
    mode  = mode
  )
  title(paste("Branch", mode, "HPDs"))
  dev.off()
}

############################################################
## 8. EXTRACT TIME SLICE FOR ANIMATION
############################################################

# Generate a sequence of time points across the tree
times <- seq(min(tr$nodeTimes), max(tr$nodeTimes), length.out = ANIMATION_FRAMES)

# Extract a single example time slice
slice <- pts_and_hpds_at_time(
  timePt = times[15],   # 15th frame as example
  tr     = tr
)

# Inspect slice: coordinates, HPDs, discrete probabilities
print(slice)

############################################################
## 9. SPATIAL ANIMATION DRIVER
############################################################

# Generates frame-by-frame PNGs for movie
# Features:
# - Interpolates node positions at each time
# - Draws HPD polygons
# - Colors by discrete trait
# - Saves frames for later ffmpeg video compilation

animate_phyloRender(
  tr,
  times   = times,
  outdir  = FRAME_DIR,
  prefix  = ANIMATION_PREFIX,
  propIndex = ANIMATION_PROP
)

message("Animation frames saved in: ", FRAME_DIR)
message("Use ffmpeg to compile: ffmpeg -r 10 -i frame_%04d.png movie.mp4")

############################################################
## 10. TEST DATA GENERATION (if no real BEAST MCC available)
############################################################

# Example: simple random tree with 20 tips
if(!file.exists(MCC_FILE)){
  set.seed(1)
  tr_test <- rtree(20)
  write.tree(tr_test, file="data/example_mcc.tree")
  message("Test MCC tree created at data/example_mcc.tree")
}

############################################################
## 11. END OF MAIN SCRIPT
############################################################

#message("MCC analysis complete. Plots and frames saved to output directories.")