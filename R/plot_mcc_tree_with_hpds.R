############################################################
## SCRIPT CREDENTIALS 
############################################################
## Tool       :  # phyloRender as part of the 'phyloAct' program
## Current script    :  # plot_mcc_tree_with_hpds.R
## Description       :  # Visualize the MCC tree in spatial layout with HPD intervals, including optional coloring by discrete traits or time.
## Author            :  # Jonathan Chan
## Affiliation       :  # Institute of Health and Community Medicine
## Contact Email     :  # jonatasp92@gmail.com
## Date Created      :  # 2024-03-16
## Last Updated      :  # 2025-12-29
## Version           :  # v1.1-beta
## R Version         :  # Recommended on R version 4.4.2
## Required Packages :  # ape, ggtree, treeio, maps, phangorn, ggplot2
## Input Files       :  # requires mcc tr in nexus format
## Output Files      :  # plots, tables, animations
#' Plot MCC tree spatially with HPDs and optional discrete or temporal colouring
plot_mcc_tree_with_hpds <- function(
  tr,                     # # The MCC tree object. Must contain $latlon and $hpds (mandatory)
  xlim=c(-180,180),       # # Longitude limits of the map. Modify to zoom in/out horizontally
  ylim=c(-60,80),         # # Latitude limits of the map. Modify to zoom in/out vertically
  tlim=c(NA,NA),          # # Time limits for nodes/edges. Use numeric range, e.g., c(0,50). NA = no restriction
  propIndex=0,            # # If >0, colors nodes by discrete trait in column `propIndex` of tr$props
  show.hpds=TRUE,         # # TRUE/FALSE: Show HPD polygons (uncertainty around node locations)
  solid.pts=TRUE,         # # TRUE/FALSE: Whether to draw solid points for nodes (currently always TRUE in code)
  show.legend=TRUE,       # # TRUE/FALSE: Show legend for discrete trait colors
  legpos="bottomleft",    # # Legend position on the plot. Can be "topright", "topleft", etc.
  useWorldHires=FALSE,    # # TRUE/FALSE: Use high-resolution world map. Slower but prettier
  new.plot=TRUE,          # # TRUE/FALSE: Create a new plot. FALSE overlays on existing plot
  bgcol="grey90",         # # Background color of the map. Change to any R color
  fcol="white",           # # Fill color for map regions
  bdcol="grey70",         # # Border color for map regions
  fill=TRUE,              # # TRUE/FALSE: Fill map regions with color
  lcex=1                  # # Legend text size. Increase/decrease to resize
) {

  stopifnot(!is.null(tr$latlon), !is.null(tr$hpds)) # # Ensure tree object has coordinates and HPDs

  # # If a new plot is requested, draw the base map
  if (new.plot) {
    if (useWorldHires) {
      maps::map("worldHires", xlim=xlim, ylim=ylim,
                col=fcol, fill=fill, border=bdcol, bg=bgcol)
    } else {
      maps::map("world", xlim=xlim, ylim=ylim,
                col=fcol, fill=fill, border=bdcol, bg=bgcol)
    }
  }

  # # Default colors for nodes and edges
  pcol  <- rep(rgb(0,0,1,0.3), nrow(tr$latlon))  # Node fill color (semi-transparent blue)
  pcol2 <- rep(rgb(0,0,1,0.6), nrow(tr$latlon))  # Node point color (darker blue)
  ecol  <- rep(rgb(0,0,0,0.4), nrow(tr$edge))    # Edge color (semi-transparent black)

  # # If propIndex > 0, recolor nodes/edges based on discrete traits
  if (propIndex > 0) {
    utraits <- tr$uprops[[propIndex]]               # Unique trait values
    cols    <- get_BEAST_cols(length(utraits))      # Colors assigned to traits
    inds    <- match(tr$props[,propIndex], utraits) # Map nodes to trait colors
    pcol    <- cols[inds]                            # Node fill colors
    pcol2   <- cols[inds]                            # Node point colors
    ecol    <- cols[match(tr$props[tr$edge[,1],propIndex], utraits)] # Edge colors

    # # Add a legend if requested
    if (show.legend) {
      graphics::legend(legpos, utraits, pch=21,
                       pt.bg=cols, bty="n", cex=lcex)
    }
  }

  # # Determine which nodes and edges to display based on tlim
  ok_nodes <- seq_len(nrow(tr$latlon))  # All nodes by default
  ok_edges <- seq_len(nrow(tr$edge))    # All edges by default

  if (!is.na(tlim[1])) {
    ok_nodes <- which(tr$nodeTimes >= tlim[1] & tr$nodeTimes <= tlim[2])
    ok_edges <- which(tr$nodeTimes[tr$edge[,2]] >= tlim[1] &
                      tr$nodeTimes[tr$edge[,1]] <= tlim[2])
  }

  # # Draw HPD polygons if requested
  if (show.hpds) {
    for (i in rev(ok_nodes)) {
      graphics::polygon(tr$hpds[[i]],
                        col=pcol[i],
                        border=NA)
    }
  }

  # # Draw arrows representing tree edges
  graphics::arrows(
    tr$latlon[tr$edge[ok_edges,1],2],
    tr$latlon[tr$edge[ok_edges,1],1],
    tr$latlon[tr$edge[ok_edges,2],2],
    tr$latlon[tr$edge[ok_edges,2],1],
    length=0.05,
    col=ecol[ok_edges]
  )

  # # Draw points for nodes
  graphics::points(tr$latlon[ok_nodes,2],
                   tr$latlon[ok_nodes,1],
                   pch=21, bg=pcol2[ok_nodes])
}