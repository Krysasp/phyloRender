############################################################
## test.R — Test Script for main.R Functions
## Requires: main.R sourced, BEAST MCC Nexus tree file
############################################################

# 0. Load main workflow functions
source("main.R")  # Make sure main.R is in the same directory

############################################################
## 1. TEST DATA INPUT
############################################################

# Path to your Nexus MCC tree file
# Replace "example_test.nexus" with your test Nexus file
TEST_TREE_FILE <- "data/example_test.nexus"

# Check if file exists
if(!file.exists(TEST_TREE_FILE)){
  stop("Please provide a Nexus MCC tree file at: ", TEST_TREE_FILE)
}

############################################################
## 2. READ MCC TREE
############################################################

# Reads MCC tree, extracts HPDs, discrete traits, and branch rates
tr <- read_beast_mcc(
  file     = TEST_TREE_FILE,
  spatial  = TRUE,   # extract lat/lon + HPDs if present
  discrete = TRUE,   # extract discrete traits if present
  rates    = TRUE    # extract branch substitution rates if present
)

# Inspect tree structure
print(tr)

############################################################
## 3. PLOT TREE TO CHECK STRUCTURE
############################################################

plot_discrete_tree(tr, propIndex = 1, legpos = "bottomleft", edgeCols = TRUE)
plot_mcc_tree_with_hpds(tr, xlim = c(-180,180), ylim = c(-60,80), propIndex = 1)

############################################################
## 4. EXTRACT TIME SLICE
############################################################

times <- seq(min(tr$nodeTimes), max(tr$nodeTimes), length.out = 5)
slice <- pts_and_hpds_at_time(timePt = times[3], tr = tr)
print(slice)

############################################################
## 5. BRANCH HPD LABELING DEMO
############################################################

plot_mcc_tree_with_branch_hpds(tr, nodes = (length(tr$tip.label)+1):(length(tr$tip.label)+3), mode="emergence")
plot_mcc_tree_with_branch_hpds(tr, nodes = (length(tr$tip.label)+1):(length(tr$tip.label)+3), mode="rate")

############################################################
## 6. ANIMATION DRIVER (optional, small test)
############################################################

animate_phyloRender(
  tr,
  times   = times,
  outdir  = "output/test_frames",
  prefix  = "test_frame",
  propIndex = 1
)

message("Test workflow complete. Check plots and animation frames in output directories.")
