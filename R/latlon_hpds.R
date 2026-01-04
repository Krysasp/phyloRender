############################################################
## FUNCTION: addLatLonHPD
############################################################
## Purpose:
##   Attach latitude and longitude HPD intervals to an MCC tree object.
##   Supports multiple HPD levels, calculates centroids, and generates
##   bounding ellipses for visualization.
##
## Parameters:
## @param tr
##   MCC tree object (ape::phylo) containing 'details' and 'latlon'.
##
## @param latlonName
##   Prefix used in metadata for latitude and longitude HPD entries.
##   Default is "latlon". Change if metadata uses custom names.
##
## @param hpdLevels
##   Integer vector of HPD levels to extract (e.g., c(1,2) for 80% and 95% HPD).
##   Default is 1 (80% HPD).
##
## @param normalizeLon
##   Logical, whether to normalize longitudes to [-180,180]. Default = TRUE.
##
## Usage Instructions:
##   - Modify latlonName if tree uses a different metadata naming convention.
##   - Add multiple HPD levels to extract multiple confidence intervals simultaneously.
##   - Generated centroids and ellipses can be used directly for plotting or animation.
############################################################
addLatLonHPD <- function(tr,
                            latlonName = "latlon",
                            hpdLevels = 1,
                            normalizeLon = TRUE) {

  # --- Initialize storage lists
  tr$latHPDs  <- vector("list", length(tr$details))
  tr$lonHPDs  <- vector("list", length(tr$details))
  tr$HPDcentroids <- matrix(NA, nrow=length(tr$details), ncol=2)
  tr$HPDellipses  <- vector("list", length(tr$details))

  # --- Loop through each node/tip
  for (i in seq_along(tr$details)) {

    # --- Storage per node for multiple HPDs
    nodeLatHPDs <- list()
    nodeLonHPDs <- list()

    for (lvl in hpdLevels) {
      # --- Regex patterns for latitude and longitude HPD
      latPattern <- paste0(latlonName,"1_80\\%HPD_", lvl, "=\\{.+\\}")
      lonPattern <- paste0(latlonName,"2_80\\%HPD_", lvl, "=\\{.+\\}")

      # --- Find matches
      latMatch <- gregexpr(latPattern, tr$details[i])[[1]]
      lonMatch <- gregexpr(lonPattern, tr$details[i])[[1]]

      # --- Extract numeric values or fallback
      latVals <- if(latMatch[1] > 0) {
        as.numeric(strsplit(gsub(".*\\{|\\}","",regmatches(tr$details[i], latMatch)), ",")[[1]])
      } else tr$latlon[i,1]

      lonVals <- if(lonMatch[1] > 0) {
        as.numeric(strsplit(gsub(".*\\{|\\}","",regmatches(tr$details[i], lonMatch)), ",")[[1]])
      } else tr$latlon[i,2]

      # --- Normalize longitude if requested
      if(normalizeLon) lonVals <- ifelse(lonVals > 180, lonVals - 360, lonVals)

      nodeLatHPDs[[lvl]] <- latVals
      nodeLonHPDs[[lvl]] <- lonVals
    }

    # --- Assign multi-HPDs
    tr$latHPDs[[i]] <- nodeLatHPDs
    tr$lonHPDs[[i]] <- nodeLonHPDs

    # --- Compute centroid of first HPD level (default for plotting/animation)
    tr$HPDcentroids[i,] <- c(mean(nodeLatHPDs[[hpdLevels[1]]], na.rm=TRUE),
                             mean(nodeLonHPDs[[hpdLevels[1]]], na.rm=TRUE))

    # --- Generate bounding ellipse for visualization (simple convex hull)
    xy <- cbind(nodeLonHPDs[[hpdLevels[1]]], nodeLatHPDs[[hpdLevels[1]]])
    if(nrow(xy) > 2) {
      ch <- chull(xy)
      tr$HPDellipses[[i]] <- xy[c(ch, ch[1]), , drop=FALSE] # closed polygon
    } else {
      tr$HPDellipses[[i]] <- xy # fallback to raw points if insufficient for hull
    }
  }

  return(tr)
}