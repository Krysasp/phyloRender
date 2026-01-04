############################################################
## FUNCTION: import_mcc_tree_with_traits
############################################################
## Description:
##   Imports a BEAST MCC tree file with [&R] annotations and extracts
##   tip and node metadata including multiple continuous traits.
##   Continuous traits (e.g., lat/lon, altitude, rates) are automatically
##   stored as matrices for downstream analysis.
##
## Features:
##   - Extracts multiple continuous traits simultaneously
##   - Automatically generates matrices for each trait in tr$trait_matrices
##   - Compatible with HPD fitting, plotting, and animation workflows
##   - Maintains tip labels, node labels, and node times for compatibility with main.R
##
## Suggested Enhancements:
##   - Filter nodes or tips based on trait thresholds
##   - Scale geographic or numeric traits for plotting or distance calculations
##
## @param treeFile
##   Path to the MCC tree (Newick/Nexus with [&R] annotations)
##
## @param traitKeys
##   Character vector of continuous trait prefixes to extract (default: "latlon")
##
## @param useDecimalDates
##   Logical flag indicating whether to convert tip dates to decimal years (default: TRUE)
##
## @param maxTipDate
##   Optional numeric specifying the youngest tip date; automatically inferred if zero
##
## @return
##   An ape::phylo object with additional fields:
##     - tr$tip.details    : original metadata strings for tips
##     - tr$node.details   : original metadata strings for nodes
##     - tr$details        : combined tip + node metadata
##     - tr$tip.label      : cleaned tip labels
##     - tr$trait_matrices : named list of matrices for each continuous trait
############################################################

read_latlon_mcc_tr <- function(treeFile,
                                        traitKeys = c("latlon"),
                                        useDecimalDates = TRUE,
                                        maxTipDate = 0) {

  # --- Read file
  treeLines <- readLines(treeFile)
  treeLine <- treeLines[grep("TREE", treeLines)]
  treeLine <- strsplit(treeLine, "\\[\\&R\\] ")[[1]][2]

  # --- Detect and extract metadata annotations
  metadataStarts <- gregexpr("\\[\\&", treeLine)[[1]]
  metadataEnds   <- gregexpr("\\]", treeLine)[[1]]
  nNodes <- length(metadataStarts)
  metadataStrings <- character(nNodes)

  for (i in seq_len(nNodes)) {
    metadataStrings[i] <- substring(treeLine, metadataStarts[i], metadataEnds[i])
    # replace with placeholder
    treeLine <- gsub(metadataStrings[i], paste0("_nodeMeta_", i), treeLine, fixed = TRUE)
  }

  # --- Read tree structure
  tree <- ape::read.tree(text = treeLine)

  # --- Map placeholders to actual metadata
  tree$tip.details  <- metadataStrings[as.integer(apply(as.matrix(tree$tip.label), 1,
                                                        getEl, ind=1, fromEnd=TRUE, sep="_"))]
  tree$node.details <- metadataStrings[as.integer(apply(as.matrix(tree$node.label), 1,
                                                        getEl, ind=1, fromEnd=TRUE, sep="_"))]

  # --- Clean tip labels
  tree$tip.label <- apply(as.matrix(tree$tip.label), 1, getEl, ind=1, sep="_")

  # --- Apply translation table if present
  translationTbl <- getTranslation(treeLines)
  tree$tip.label <- translationTbl[match(tree$tip.label, translationTbl[,1]),2]

  # --- Remove node labels for compatibility
  tree$node.label <- NULL

  # --- Combine all node metadata
  tree$details <- c(tree$tip.details, tree$node.details)

  # --- Compute decimal tip dates if requested
  if (useDecimalDates) {
    tipDates <- as.numeric(apply(as.matrix(tree$tip.label), 1,
                                 getEl, ind=1, sep="\\|", fromEnd=TRUE))
    maxTipDate <- max(tipDates)
  }

  # --- Compute node times
  tree <- nodeTimes(tree, youngestTip = maxTipDate)

  # --- Extract continuous traits
  tree$trait_matrices <- list()
  for (trait in traitKeys) {

    # Attempt to extract two components (e.g., lat/lon)
    comp1 <- apply(as.matrix(tree$details), 1, getEl, ind=2, sep=paste0(trait,"1="))
    comp2 <- apply(as.matrix(tree$details), 1, getEl, ind=2, sep=paste0(trait,"2="))

    comp1 <- suppressWarnings(as.numeric(apply(as.matrix(comp1), 1, getEl, ind=1, sep=",")))
    comp2 <- suppressWarnings(as.numeric(apply(as.matrix(comp2), 1, getEl, ind=1, sep=",")))

    # Generate matrix
    if (all(is.na(comp2))) {
      # Single-column trait
      tree$trait_matrices[[trait]] <- matrix(comp1, ncol = 1)
      colnames(tree$trait_matrices[[trait]]) <- trait
    } else {
      # Two-column trait
      tree$trait_matrices[[trait]] <- cbind(comp1, comp2)
      colnames(tree$trait_matrices[[trait]]) <- paste0(trait, c("_1","_2"))
    }
  }

  return(tree)
}
