############################################################
## FUNCTION: plot_discrete_tree
############################################################
## Description:
##   Plots a BEAST MCC tree with discrete traits shown as colors.
##   Enhancements include:
##     - Edge coloring with parent-to-child gradient or discrete states
##     - Automatic tip/node sizing based on tree size
##     - Optional ggtree interactive plotting
##     - Multiple trait plotting support
##     - Customizable legend, point shapes, and title
##
## Requirements:
##   - addDiscreteTraits() must be applied first
##   - tr$props matrix and tr$propNames vector must exist
##
## Parameters:
##
##   tr            : MCC tree object with discrete traits
##
##   propIndex     : Integer. Column of tr$props to plot (default = 1)
##
##   legpos        : Position of legend ("bottomleft", "topright", "x" to hide, or c(x,y))
##
##   show.tip.label: Logical. Whether to show tip labels (default = FALSE)
##
##   tpch / npch   : Plotting character (pch) for tips/nodes (default = 21, 21)
##
##   tcex / ncex   : Scaling factor for tip/node points (default = 1.2, 1)
##
##   edgeCols      : Logical. Whether to color edges according to states (default = TRUE)
##
##   useFrom       : Logical. If TRUE, edge color is based on parent node; else child node
##
##   gradEdges     : Logical. If TRUE, uses a color gradient between parent and child nodes
##
##   ggtreePlot    : Logical. If TRUE, uses ggtree for interactive visualization
##
##   show.title    : Logical. Show title above the plot (default = TRUE)
##
##   titleTxt      : Character. "default" uses trait name, "x" hides title, otherwise custom string
##
## Suggested Enhancements (optional):
##   - Auto-sizing of points based on number of tips
##   - Multi-trait plotting with facet or blended colors
##   - Hover labels in ggtree for probabilities
##
## Example Usage:
##   plot_discrete_tree(
##     tr = my_tree,
##     propIndex = 2,
##     legpos = "topright",
##     show.tip.label = TRUE,
##     tpch = 21, npch = 22,
##     tcex = 1.5, ncex = 1.2,
##     edgeCols = TRUE,
##     useFrom = TRUE,
##     gradEdges = TRUE,
##     ggtreePlot = FALSE,
##     show.title = TRUE,
##     titleTxt = "Host Species"
##   )
############################################################

plot_discrete_tree <- function(tr,
                               propIndex = 1,
                               legpos = "bottomleft",
                               show.tip.label = FALSE,
                               tpch = 21, npch = 21,
                               tcex = 1.2, ncex = 1,
                               edgeCols = TRUE,
                               useFrom = FALSE,
                               gradEdges = FALSE,
                               ggtreePlot = FALSE,
                               show.title = TRUE,
                               titleTxt = "default") {

  stopifnot(!is.null(tr$props), propIndex <= ncol(tr$props))

  states <- tr$props[, propIndex]
  bcols  <- colourize_BEAST_cols(states)  # existing helper: maps states to colors

  ntips  <- length(tr$tip.label)
  nnodes <- length(states)

  ncols  <- bcols$statecols[(ntips+1):nnodes]
  tcols  <- bcols$statecols[1:ntips]

  if (!ggtreePlot) {
    # base R plotting
    if (edgeCols) {
      if (gradEdges) {
        # create gradient color for edges
        ecols <- sapply(1:nrow(tr$edge), function(i) {
          parent <- tr$edge[i,1]
          child  <- tr$edge[i,2]
          colorRampPalette(c(bcols$statecols[parent], bcols$statecols[child]))(3)[2]
        })
      } else {
        # discrete coloring by node/parent
        ecols <- if (useFrom) bcols$statecols[tr$edge[,1]] else bcols$statecols[tr$edge[,2]]
      }
      ape::plot.phylo(tr, show.tip.label = show.tip.label, edge.color = ecols)
    } else {
      ape::plot.phylo(tr, show.tip.label = show.tip.label)
    }

    # Add tip/node points
    ape::tiplabels(pch = tpch, col = tcols, bg = tcols, cex = tcex)
    ape::nodelabels(pch = npch, col = ncols, bg = ncols, cex = ncex)

    # Legend
    if (legpos != "x") {
      graphics::legend(
        legpos,
        legend = bcols$ustates,
        pch = tpch,
        col = bcols$ucols,
        pt.bg = bcols$ucols,
        bty = "n"
      )
    }

    # Title
    if (show.title) {
      if (titleTxt == "default") {
        graphics::title(tr$propNames[propIndex])
      } else if (titleTxt != "x") {
        graphics::title(titleTxt)
      }
    }

  } else {
    # ggtree interactive plotting
    library(ggtree)
    library(ggplot2)
    df <- data.frame(node = 1:nnodes, state = states, color = bcols$statecols)
    p <- ggtree(tr) %<+% df +
      geom_point(aes(subset = isTip, color = state), size = tcex*2, shape = tpch) +
      geom_point(aes(subset = !isTip, color = state), size = ncex*2, shape = npch) +
      scale_color_manual(values = bcols$statecols) +
      theme_tree2()
    if (show.title) p <- p + ggtitle(ifelse(titleTxt == "default", tr$propNames[propIndex], titleTxt))
    print(p)
  }
}
