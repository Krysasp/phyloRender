############################################################
## SCRIPT CREDENTIALS 
############################################################
## Tool       :  # phyloRender as part of the 'phyloAct' program
## Current script    :  # read_mcc_tr.R
## Author            :  # Jonathan Chan
## Affiliation       :  # Institute of Health and Community Medicine
## Last Updated      :  # 2025-12-16
############################################################
## FUNCTION: mcc_import_and_read
############################################################
## Purpose:
##   Import a BEAST MCC tree (BEAST1 or BEAST2), extract metadata,
##   compute node times, add posterior/support values, and optionally:
##     - Highlight branches by chosen metadata column
##     - Color tip points based on metadata
##     - Show/hide tip points
##     - Switch tree layout (rectangular/cladogram/radial/unrooted)
##     - Export node and tip metadata with support to CSV
##
## Parameters:
## @param treeFile
##   Path to the MCC tree file (Newick/Nexus)
##
## @param dateSep
##   Separator in tip labels to parse decimal dates (default: "\\|")
##
## @param isBEAST2
##   TRUE if tree is from BEAST2, FALSE for BEAST1
##
## @param highlightColumn
##   Name of metadata column used to highlight branches (optional)
##
## @param tipColorColumn
##   Name of metadata column used to color tip points (optional)
##
## @param showTips
##   TRUE/FALSE to show tip points in plots
##
## @param treeStyle
##   One of "rectangular", "cladogram", "radial", "unrooted" (default: "rectangular")
##
## @param exportCSV
##   File path to export node & tip metadata with posterior/support
##   If NULL, CSV is not exported
##
## Returns:
##   A list containing:
##     - tr: ape::phylo object with additional fields
##     - tipColors: vector of colors for tips (if tipColorColumn used)
##     - highlightedEdges: vector of edge indices that match highlightColumn
############################################################

read_mcc_tr <- function(treeFile,
                                dateSep = "\\|",
                                isBEAST2 = FALSE,
                                highlightColumn = NULL,
                                tipColorColumn = NULL,
                                showTips = TRUE,
                                treeStyle = "rectangular",
                                exportCSV = NULL) {

  # --- Read raw tree
  treeLines <- readLines(treeFile)
  rawTree <- treeLines[length(treeLines)-1]
  treeStr <- if (isBEAST2) {
    strsplit(rawTree, "tree TREE1 = ")[[1]][2]
  } else {
    strsplit(rawTree, "\\[\\&R\\] ")[[1]][2]
  }

  # --- Extract node metadata
  metaStart <- gregexpr("\\[\\&", treeStr)[[1]]
  metaEnd   <- gregexpr("\\]", treeStr)[[1]]
  nNodes <- length(metaStart)
  nodeMeta <- character(nNodes)

  for (i in seq_len(nNodes)) {
    nodeMeta[i] <- substring(treeStr, metaStart[i], metaEnd[i])
    treeStr <- gsub(nodeMeta[i], paste0("_nodeMeta_", i), treeStr, fixed = TRUE)
  }

  # --- Read tree structure
  tr <- ape::read.tree(text = treeStr)

  # --- Map metadata to tips and nodes
  tr$tip.details  <- nodeMeta[as.integer(apply(as.matrix(tr$tip.label), 1,
                                               getEl, ind=1, fromEnd=TRUE, sep="_"))]
  tr$node.details <- nodeMeta[as.integer(apply(as.matrix(tr$node.label), 1,
                                               getEl, ind=1, fromEnd=TRUE, sep="_"))]

  # --- Clean tip labels
  tr$tip.label <- apply(as.matrix(tr$tip.label), 1, getEl, ind=1, sep="_")

  # --- Translation table
  translationTbl <- if (isBEAST2) getTranslation_BEAST2(treeLines)
                    else getTranslation(treeLines)
  tr$tip.label <- translationTbl[match(tr$tip.label, translationTbl[,1]),2]

  # --- Remove node labels & combine metadata
  tr$node.label <- NULL
  tr$details <- c(tr$tip.details, tr$node.details)

  # --- Compute tip dates and node times
  tipDates <- as.numeric(apply(as.matrix(tr$tip.label), 1,
                               getEl, ind=1, sep=dateSep, fromEnd=TRUE))
  youngestTip <- max(tipDates)
  tr <- nodeTimes(tr, youngestTip = youngestTip)

  # --- Add posterior probabilities and heights
  tr <- addPosterior(tr)
  if (!isBEAST2) tr <- addHeights(tr)

  # --- Prepare tip coloring
  tipColors <- NULL
  if (!is.null(tipColorColumn)) {
    colVals <- sapply(tr$tip.details, function(md) {
      getEl(md, ind=tipColorColumn, sep=",")
    })
    tipColors <- rainbow(length(unique(colVals)))[as.numeric(factor(colVals))]
  }

  # --- Highlight edges based on metadata column
  highlightedEdges <- NULL
  if (!is.null(highlightColumn)) {
    highlightedEdges <- which(sapply(seq_len(nrow(tr$edge)), function(i) {
      val <- getEl(tr$node.details[tr$edge[i,2]], ind=highlightColumn, sep=",")
      !is.na(val) && nchar(val) > 0
    }))
  }

  # --- Tree layout for plotting
  layoutType <- switch(tolower(treeStyle),
                       rectangular = "phylogram",
                       cladogram   = "cladogram",
                       radial      = "fan",
                       unrooted    = "unrooted",
                       "phylogram") # default

  # --- Export metadata & support
  if (!is.null(exportCSV)) {
    nodeSupport <- sapply(seq_len(tr$Nnode), function(i) tr$posterior[i])
    df <- data.frame(
      node = (length(tr$tip.label)+1):(length(tr$tip.label)+tr$Nnode),
      metadata = tr$node.details,
      posterior = nodeSupport
    )
    write.csv(df, file=exportCSV, row.names=FALSE)
  }

  return(list(
    tr = tr,
    tipColors = tipColors,
    highlightedEdges = highlightedEdges,
    layout = layoutType,
    showTips = showTips
  ))
}
