############################################################
## FUNCTION: swapLatLonAdvanced
############################################################
## Purpose:
##   Flexibly flip/interchange latitude and longitude
##   in MCC tree objects. Supports tip coordinates, HPD polygons,
##   and multiple HPD levels. Optional visualization for sanity checks.
##
## Parameters:
## @param tr
##   MCC tree object containing:
##     - tr$latlon        : tip lat/lon matrix
##     - tr$lat80 / tr$lon80 : list of HPD lat/lon per node
##     - tr$hpds           : node HPD polygons/points
##     - tr$hpd_diffs      : edge HPD differences
##
## @param swapTips (logical, default TRUE)
##   If TRUE, swaps lat/lon of tip coordinates
##
## @param swapHPDs (logical, default TRUE)
##   If TRUE, swaps lat/lon of node HPD points/polygons
##
## @param hpdLevels (integer vector or NULL, default NULL)
##   Which HPD levels to swap (applies if multiple levels exist)
##   NULL swaps all available HPD lists.
##
## @param plotCheck (logical, default FALSE)
##   If TRUE, plots before/after swap for sanity check
##
## Usage:
##   tr <- switch_latlon_advanced(tr, swapTips=TRUE, swapHPDs=TRUE,
##                                hpdLevels=NULL, plotCheck=FALSE)
##
## Suggested Enhancements:
##   - Add option to swap only specific nodes or tips by index.
##   - Normalize longitude values after swapping if needed.
##   - Can integrate with animation/plot scripts to auto-flip coordinates.
############################################################
# --- Helper function to swap columns of a 2-column matrix
switch_list_el <- function(mat) {
  if(!is.matrix(mat) || ncol(mat) != 2) stop("Input must be a 2-column matrix")
  mat[, c(2,1)]
}

# --- Main advanced swapping function
switch_latlon_advanced <- function(tr,
                                   swapTips=TRUE,
                                   swapHPDs=TRUE,
                                   hpdLevels=NULL,
                                   plotCheck=FALSE) {

  # --- Optional pre-swap plot
  if(plotCheck && swapTips && !is.null(tr$latlon)) {
    plot(tr$latlon[,1], tr$latlon[,2], main="Tip coordinates before swap",
         xlab="Lat", ylab="Lon", pch=21, bg="red")
  }

  # --- Swap tip coordinates
  if(swapTips && !is.null(tr$latlon)) {
    tr$latlon <- switch_list_elements(tr$latlon)
  }

  # --- Swap HPD latitude/longitude lists
  if(swapHPDs) {
    # Decide which HPD levels to swap
    if(is.null(hpdLevels)) hpdLevels <- seq_along(tr$lat80)
    
    for(h in hpdLevels) {
      if(!is.null(tr$lat80[[h]]) & !is.null(tr$lon80[[h]])) {
        tmp        <- tr$lat80[[h]]
        tr$lat80[[h]] <- tr$lon80[[h]]
        tr$lon80[[h]] <- tmp
      }
    }

    # Swap fitted HPD polygons and edge differences
    if(!is.null(tr$hpds)) {
      tr$hpds <- lapply(tr$hpds, switch_list_elements)
    }
    if(!is.null(tr$hpd_diffs)) {
      tr$hpd_diffs <- lapply(tr$hpd_diffs, switch_list_elements)
    }
  }

  # --- Optional post-swap plot
  if(plotCheck && swapTips && !is.null(tr$latlon)) {
    plot(tr$latlon[,1], tr$latlon[,2], main="Tip coordinates after swap",
         xlab="Lat", ylab="Lon", pch=21, bg="blue")
  }

  # --- Return modified tree
  return(tr)
}
