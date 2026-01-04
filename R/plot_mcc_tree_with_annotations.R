############################################################
## SCRIPT CREDENTIALS 
############################################################
## Tool       :  # phyloRender as part of the 'phyloAct' program
## Current script    :  # plot_mcc_tree_with_annotations.R
## Author            :  # Jonathan Chan
## Affiliation       :  # Institute of Health and Community Medicine
## Last Updated      :  # 2025-11-28
############################################################
## FUNCTION: plot_mcc_tree_with_annotations
############################################################
## Description:
##   Enhanced MCC tree plotting function with branch-level HPD and rate annotations.
##   Features:
##     - Emergence time HPDs (decimal years)
##     - Substitution rate HPDs
##     - Automatic node selection based on HPD width or rate thresholds
##     - Color-coded labels for easy visualization
##     - Optional simultaneous plotting of both emergence and rate HPDs
##     - Compatible with ape::plot.phylo and optionally ggtree
##
## Requirements:
##   - tr$nodeTimes must exist
##   - For 'emergence': tr$height.lower95 and tr$height.upper95
##   - For 'rate': tr$rate_range_1 and tr$rate_range_2
##
## Parameters:
##
##   tr            : MCC tree object with nodeTimes and optional HPDs/rates
##
##   nodes         : Vector of node numbers to annotate. 
##                   Default: auto-selected nodes exceeding thresholds
##
##   mode          : Character vector. Either "emergence", "rate", or "both"
##
##   digits        : Number of decimal places for label rounding (default = 3)
##
##   offset        : Horizontal offset for label placement (default = 0.02)
##
##   hpd_thresh    : Numeric. Minimum HPD width to trigger auto-node selection (default = 0.5)
##
##   rate_thresh   : Numeric. Minimum rate to trigger auto-node selection (default = 1)
##
##   col_high      : Color for high HPD width or high rates (default = "red")
##
##   col_low       : Color for low HPD width or low rates (default = "black")
##
## Suggested Enhancements:
##   - Color branches according to rate/HPD values
##   - Combine with plot_discrete_tree() for simultaneous discrete trait visualization
##   - ggtree interactive plotting for hoverable labels
##
## Usage Example:
##   plot_mcc_tree_with_annotations(
##       tr = my_tree,
##       nodes = NULL,          # auto-select nodes
##       mode = "both",
##       digits = 2,
##       offset = 0.05,
##       hpd_thresh = 0.3,
##       rate_thresh = 1.2
##   )
############################################################

plot_mcc_tree_with_annotations <- function(
  tr,
  nodes = NULL,
  mode = c("emergence","rate","both"),
  digits = 3,
  offset = 0.02,
  hpd_thresh = 0.5,
  rate_thresh = 1,
  col_high = "red",
  col_low = "black"
) {
  
  mode <- match.arg(mode)
  
  stopifnot(!is.null(tr$nodeTimes))
  if(mode %in% c("emergence","both")) stopifnot(!is.null(tr$height.lower95), !is.null(tr$height.upper95))
  if(mode %in% c("rate","both")) stopifnot(!is.null(tr$rate_range_1), !is.null(tr$rate_range_2))
  
  # Automatic node selection if nodes not provided
  if(is.null(nodes)) {
    nodes <- seq_along(tr$nodeTimes)
    if(mode %in% c("emergence","both")) {
      nodes <- nodes[ (tr$height.upper95 - tr$height.lower95) > hpd_thresh ]
    }
    if(mode %in% c("rate","both")) {
      nodes <- unique(c(nodes, which(tr$rate_range_2 > rate_thresh)))
    }
  }
  
  # Plot tree
  ape::plot.phylo(tr, show.tip.label = FALSE)
  
  # Loop through nodes for annotation
  for(n in nodes) {
    parent <- tr$edge[which(tr$edge[,2]==n),1]
    x <- mean(tr$nodeTimes[c(parent,n)])
    y <- n
    
    # Build label depending on mode
    lab <- ""
    col <- col_low
    if(mode == "emergence") {
      lab <- paste0(
        round(tr$height.lower95[n],digits), "–",
        round(tr$height.upper95[n],digits)
      )
      col <- if((tr$height.upper95[n]-tr$height.lower95[n]) > hpd_thresh) col_high else col_low
    } else if(mode == "rate") {
      lab <- paste0(
        signif(tr$rate_range_1[n],digits), "–",
        signif(tr$rate_range_2[n],digits)
      )
      col <- if(tr$rate_range_2[n] > rate_thresh) col_high else col_low
    } else if(mode == "both") {
      lab <- paste0(
        round(tr$height.lower95[n],digits), "–", round(tr$height.upper95[n],digits),
        " / ",
        signif(tr$rate_range_1[n],digits), "–", signif(tr$rate_range_2[n],digits)
      )
      col <- if((tr$height.upper95[n]-tr$height.lower95[n]) > hpd_thresh | tr$rate_range_2[n] > rate_thresh) col_high else col_low
    }
    
    graphics::text(x + offset, y, labels = lab, col = col, cex = 0.7)
  }
  
  # Optional: ggtree interactive plotting (requires ggtree)
  # library(ggtree)
  # ggtree::ggtree(tr) + ggtree::geom_text(aes(label=lab, x=x+offset, y=y, color=col))
}
