############################################################
## SCRIPT CREDENTIALS 
############################################################
## Tool       :  # phyloRender as part of the 'phyloAct' program
## Current script    :  # hpds_fit_standard.R
## Author            :  # Jonathan Chan
## Affiliation       :  # Institute of Health and Community Medicine
## Last Updated      :  # 2025-12-17
############################################################
## FUNCTION: fit_HPDs_to_standard
############################################################
## Description:
##   Standardizes HPD intervals (lat/lon) for tree branches/nodes.
##   Can handle multiple HPD levels (80%, 95%) and automatically
##   fits ellipses to spatial distributions using 'conicfit'.
##
## Features:
##   - Supports multiple HPD levels: tr$lat80/lon80, tr$lat95/lon95, etc.
##   - Computes branch-level differences for animation (tr$hpd_diffs)
##   - Handles small ranges and missing data robustly
##   - Optional parameters for number of points and minimal tolerance
##
## Suggested Enhancements:
##   - Automatically scale coordinates to km
##   - Combine with continuous traits for animated uncertainty visualization
##   - Use for plotting or spatial phylodynamics animations
##
## @param tr
##   MCC tree object containing lat/lon HPD lists for each level.
##
## @param levels
##   Character vector of HPD levels to fit (default: c("80","95"))
##   Corresponding tree fields should exist: e.g., lat80/lon80, lat95/lon95
##
## @param npts
##   Number of points to generate for ellipse/polygon (default: 50)
##
## @param ltol
##   Minimum difference threshold to prevent zero-size HPDs (default: 0.005)
##
## @return
##   Updated tree object with fields for each HPD level:
##     - tr$hpds_LEVEL   : fitted HPD polygons for each node
##     - tr$hpd_diffs_LEVEL : branch differences (to - from)
##
## Usage:
##   tr <- fit_HPDs_to_standard(tr, levels=c("80","95"), npts=100)
############################################################

fit_HPDs_to_standard <- function(tr, levels=c("80","95"), npts=50, ltol=0.005) {

  # Loop through each specified HPD level
  for (level in levels) {
    lat_field <- paste0("lat", level)
    lon_field <- paste0("lon", level)

    if (!(lat_field %in% names(tr)) | !(lon_field %in% names(tr))) {
      warning(paste0("HPD fields ", lat_field, "/", lon_field, " not found. Skipping."))
      next
    }

    lat_list <- tr[[lat_field]]
    lon_list <- tr[[lon_field]]
    n_nodes <- length(lat_list)

    # Prepare lists to store fitted polygons and branch differences
    hpds <- vector("list", n_nodes)
    hpd_diffs <- vector("list", length(tr$edge.length))

    # Fit HPD polygons for each node/branch
    for (i in 1:n_nodes) {
      ytemp <- lat_list[[i]]
      xtemp <- lon_list[[i]]

      # Check for sufficient variation
      xgood <- length(xtemp) > 3 && (max(xtemp, na.rm=TRUE)-min(xtemp, na.rm=TRUE)) > ltol
      ygood <- length(ytemp) > 3 && (max(ytemp, na.rm=TRUE)-min(ytemp, na.rm=TRUE)) > ltol

      if (xgood & ygood) {
        kk <- which(is.finite(xtemp) & is.finite(ytemp))
        xy <- cbind(xtemp[kk], ytemp[kk])

        # Remove duplicate last points
        nk <- nrow(xy)
        if (nk > 3 && (all(xy[nk,] == xy[nk-1,]) | all(xy[nk,] == xy[1,]))) {
          xy <- xy[1:(nk-1),]
        }

        ellipDirect <- conicfit::EllipseDirectFit(xy)
        if (all(is.na(ellipDirect))) {
          ellipG <- c(mean(xtemp, na.rm=TRUE), mean(ytemp, na.rm=TRUE), ltol, ltol, 0)
        } else {
          ellipG <- conicfit::AtoG(ellipDirect)$ParG
        }

      } else {
        # Handle small ranges
        xmid <- mean(xtemp, na.rm=TRUE)
        ymid <- mean(ytemp, na.rm=TRUE)
        xrange <- if (xgood) (max(xtemp, na.rm=TRUE)-min(xtemp, na.rm=TRUE))/2 else ltol
        yrange <- if (ygood) (max(ytemp, na.rm=TRUE)-min(ytemp, na.rm=TRUE))/2 else ltol
        ellipG <- c(xmid, ymid, xrange, yrange, 0)
      }

      # Generate fitted polygon points
      hpds[[i]] <- conicfit::calculateEllipse(ellipG[1], ellipG[2], ellipG[3], ellipG[4], 180/pi*ellipG[5], steps=npts)
    }

    # Compute branch-level differences
    for (j in 1:length(tr$edge.length)) {
      fromHPD <- hpds[[ tr$edge[j,1] ]]
      toHPD <- hpds[[ tr$edge[j,2] ]]
      hpd_diffs[[j]] <- toHPD - fromHPD
    }

    # Save results in tree object
    tr[[paste0("hpds_", level)]] <- hpds
    tr[[paste0("hpd_diffs_", level)]] <- hpd_diffs
  }

  return(tr)
}
