############################################################
## SCRIPT CREDENTIALS 
############################################################
## Tool       :  # phyloRender as part of the 'phyloAct' program
## Current script    :  # read_beast_mcc.R
## Author            :  # Jonathan Chan
## Affiliation       :  # Institute of Health and Community Medicine
## Last Updated      :  # 2025-12-19
############################################################
## FUNCTION: plot_mcc_tree_hpds_enhanced
############################################################
## Description:
##   Plots an MCC tree and annotates branches with HPD intervals.
##   Features include:
##     - Automatic node selection based on HPD width or rate thresholds
##     - Simultaneous plotting of emergence and rate HPDs
##     - Color-coded branches based on HPD width or rate
##     - Optional legends for clarity
##     - Optional interactive plotting with ggtree
##
## Requirements:
##   - tr: must be processed with addRates() and addHeights()
##   - ape, ggtree, and RColorBrewer packages
##
## Parameters:
##   tr            : Phylo object (MCC tree with metadata)
##   mode          : Vector of "emergence" and/or "rate" to indicate which HPDs to plot
##   auto_nodes    : Logical. If TRUE, nodes are auto-selected based on thresholds
##   hpd_width_cut : Numeric. Minimum HPD width (emergence) to annotate if auto_nodes = TRUE
##   rate_cut      : Numeric. Minimum substitution rate (or rate difference) to annotate
##   digits        : Number of digits for text labels
##   offset        : Horizontal offset for text labels
##   cex           : Text size for labels
##   use_ggtree    : Logical. If TRUE, use ggtree for interactive visualization
##   col_palette   : Vector of colors to map HPD width or rate
##
## Return:
##   - Plot of the MCC tree with annotated branches
##
## Usage Example:
##   plot_mcc_tree_hpds_enhanced(
##     tr = my_tree,
##     mode = c("emergence","rate"),
##     auto_nodes = TRUE,
##     hpd_width_cut = 0.5,
##     rate_cut = 0.2,
##     digits = 2,
##     offset = 0.03,
##     cex = 0.8,
##     use_ggtree = TRUE
##   )
##
## Notes:
##   - If auto_nodes = FALSE, user can manually specify nodes via 'nodes' vector
##   - Colors are mapped to HPD width for emergence or to mean rate for rate HPDs
##   - Requires tr to have nodeTimes, height.lower95/upper95, rate_range_1/2
############################################################

plot_mcc_tree_hpds_enhanced <- function(
  tr,
  mode = c("emergence", "rate"),
  auto_nodes = TRUE,
  hpd_width_cut = 0.5,
  rate_cut = 0.2,
  digits = 3,
  offset = 0.02,
  cex = 0.7,
  use_ggtree = FALSE,
  col_palette = RColorBrewer::brewer.pal(9, "YlOrRd")
) {

  mode <- match.arg(mode, several.ok = TRUE)

  # Check required metadata
  stopifnot(!is.null(tr$nodeTimes))
  if ("emergence" %in% mode)
    stopifnot(!is.null(tr$height.lower95), !is.null(tr$height.upper95))
  if ("rate" %in% mode)
    stopifnot(!is.null(tr$rate_range_1), !is.null(tr$rate_range_2))

  # Automatically select nodes if requested
  nodes_to_label <- NULL
  if (auto_nodes) {
    if ("emergence" %in% mode) {
      width <- tr$height.upper95 - tr$height.lower95
      nodes_to_label <- which(width >= hpd_width_cut)
    }
    if ("rate" %in% mode) {
      rate_diff <- tr$rate_range_2 - tr$rate_range_1
      rate_nodes <- which(rate_diff >= rate_cut)
      nodes_to_label <- sort(unique(c(nodes_to_label, rate_nodes)))
    }
  } else {
    nodes_to_label <- 1:Nnode(tr) # default: label all nodes
  }

  # Assign colors for branches
  node_colors <- rep("black", length(tr$edge[,2]))
  if ("emergence" %in% mode) {
    width <- tr$height.upper95 - tr$height.lower95
    cols <- grDevices::colorRampPalette(col_palette)(100)
    node_colors <- cols[as.numeric(cut(width, breaks = 100))]
  }
  if ("rate" %in% mode) {
    rate_mean <- (tr$rate_range_1 + tr$rate_range_2)/2
    cols <- grDevices::colorRampPalette(col_palette)(100)
    node_colors <- cols[as.numeric(cut(rate_mean, breaks = 100))]
  }

  # Plotting
  if (!use_ggtree) {
    # Base ape plot
    ape::plot.phylo(tr, show.tip.label = FALSE, edge.color = node_colors)
    for (n in nodes_to_label) {
      ed <- which(tr$edge[,2] == n)
      parent <- tr$edge[ed,1]
      x <- mean(tr$nodeTimes[c(parent,n)])
      y <- n
      lab <- c()
      if ("emergence" %in% mode) {
        lab <- c(lab,
          paste0("E:", round(tr$height.lower95[n],digits),
                 "–", round(tr$height.upper95[n],digits))
        )
      }
      if ("rate" %in% mode) {
        lab <- c(lab,
          paste0("R:", signif(tr$rate_range_1[n],digits),
                 "–", signif(tr$rate_range_2[n],digits))
        )
      }
      graphics::text(x + offset, y, paste(lab, collapse="\n"), cex=cex)
    }
    # Legend
    legend("topright", legend = c("Emergence HPD","Rate HPD"),
           col = c("red","blue"), pch = 15, bty = "n")
  } else {
    # ggtree interactive plot
    library(ggtree)
    g <- ggtree::ggtree(tr)
    if ("emergence" %in% mode) {
      width <- tr$height.upper95 - tr$height.lower95
      g <- ggtree::geom_tree(g, color = scales::col_numeric(col_palette, domain = width)(width))
    }
    ggtree::geom_text2(g, aes(subset = (node %in% nodes_to_label),
                              label = paste0(
                                ifelse("emergence" %in% mode,
                                       paste0("E:", round(tr$height.lower95, digits),
                                              "–", round(tr$height.upper95, digits)), ""),
                                ifelse("rate" %in% mode,
                                       paste0("\nR:", signif(tr$rate_range_1, digits),
                                              "–", signif(tr$rate_range_2, digits)), "")
                              )
    ), hjust = -0.1, size = cex*3)
  }

  invisible(nodes_to_label)
}
