############################################################
## SCRIPT CREDENTIALS 
############################################################
## Tool       :  # phyloRender as part of the 'phyloAct' program
## Current script    :  # read_beast_mcc.R
## Author            :  # Jonathan Chan
## Affiliation       :  # Institute of Health and Community Medicine
## Last Updated      :  # 2025-11-29
############################################################
## FUNCTION: read_beast_mcc
############################################################
## Description:
##   Unified and production-ready reader for BEAST MCC trees.
##   This function imports a MCC tree (Nexus or BEAST output) and
##   annotates it with spatial coordinates, discrete traits, branch rates,
##   and optionally prepares it for animation or export.
##
## Parameters:
##   file          : Character string. Path to your MCC tree file.
##                   Must be readable by read_mcc_treeio().
##
##   spatial       : Logical (TRUE/FALSE). Default TRUE.
##                   Adds latitude/longitude HPDs to tree tips and
##                   fits HPDs to a standard format for visualization.
##
##   discrete      : Logical (TRUE/FALSE). Default TRUE.
##                   Adds discrete trait annotations to tips or nodes.
##
##   rates         : Logical (TRUE/FALSE). Default TRUE.
##                   Calculates branch-specific rates and labels
##                   emergence HPDs.
##
##   animate       : Logical (TRUE/FALSE). Default FALSE.
##                   If TRUE, prepares the tree with time-sliced frames
##                   for gganimate visualization.
##
##   export_path   : Character string or NULL. Default NULL.
##                   If provided, exports annotated tree to CSV (HPDs/rates)
##                   and Newick formats.
##
##   color_by      : Character string. Default NULL.
##                   Column name (trait or rate) to color tips/branches
##                   in plotting.
##
## Usage Example:
##   tr <- read_beast_mcc(
##            file = "data/example_mcc.nexus",
##            spatial = TRUE,
##            discrete = TRUE,
##            rates = TRUE,
##            animate = TRUE,
##            export_path = "results/example_mcc",
##            color_by = "host_type"
##        )
##
## Suggested Enhancements:
##   - Filter nodes by posterior probability or HPD coverage
##   - Support reading multiple tree files into a merged object
##   - Integrate with interactive ggtree Shiny apps
##
## Notes:
##   - Requires packages: treeio, ape, ggtree, dplyr
##   - Custom functions: addLatLonHPD, fit_HPDs_to_standard,
##     addDiscreteTraits, addRates, prepare_animation_frames
############################################################

read_beast_mcc <- function(file,
                           spatial = TRUE,
                           discrete = TRUE,
                           rates = TRUE,
                           animate = FALSE,
                           export_path = NULL,
                           color_by = NULL) {
  
  # 1. Read MCC tree
  tr <- read_mcc_treeio(file)
  
  # 2. Add spatial HPDs if requested
  if (spatial) {
    tr <- addLatLonHPD(tr)
    tr <- fit_HPDs_to_standard(tr)
  }
  
  # 3. Add discrete traits
  if (discrete) {
    tr <- addDiscreteTraits(tr)
  }
  
  # 4. Add branch rates / emergence HPDs
  if (rates) {
    tr <- addRates(tr)
  }
  
  # 5. Prepare for animation if requested
  if (animate) {
    tr <- prepare_animation_frames(tr)
  }
  
  # 6. Optional coloring based on trait/rate
  if (!is.null(color_by)) {
    if (!(color_by %in% colnames(tr@data))) {
      warning(paste("color_by column", color_by, "not found in tree metadata"))
    } else {
      tr$color <- tr@data[[color_by]]
    }
  }
  
  # 7. Optional export
  if (!is.null(export_path)) {
    dir.create(dirname(export_path), showWarnings = FALSE, recursive = TRUE)
    
    # Export Newick
    write.tree(as.phylo(tr), file = paste0(export_path, ".nwk"))
    
    # Export CSV with tip data
    tip_data <- tr@data
    write.csv(tip_data, file = paste0(export_path, "_tips.csv"), row.names = FALSE)
  }
  
  return(tr)
}

############################################################
## Required supporting functions (must be defined elsewhere)
## - addLatLonHPD(tree): attach latitude/longitude HPDs to tips
## - fit_HPDs_to_standard(tree): normalize HPDs for plotting
## - addDiscreteTraits(tree): add discrete trait metadata
## - addRates(tree): compute branch rates and emergence HPDs
## - prepare_animation_frames(tree): expand tree into time slices
############################################################