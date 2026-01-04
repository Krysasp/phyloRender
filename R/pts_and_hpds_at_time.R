############################################################
## Tool       :  # phyloRender as part of the 'phyloAct' program
## Current script    :  # pts_and_hpds_at_time.R
## Author            :  # Jonathan Chan
## Affiliation       :  # Institute of Health and Community Medicine
## Last Updated      :  # 2025-11-28
############################################################
## FUNCTION: pts_and_hpds_at_time
############################################################
## Description:
##   Extracts the spatial positions, HPD polygons, and optional discrete traits
##   for a MCC tree at a specified time point.
##   Supports interpolation along branches crossing `timePt`.
##   Useful for animations, tree slices, or spatial projections.
##
## Features:
##   - Auto-selects branches intersecting `timePt`
##   - Optionally excludes specified tips
##   - Interpolates HPD polygons along branches
##   - Interpolates discrete traits if tr$fullset exists
##   - Supports color-coding HPDs by width or trait
##   - Returns both lists and data frames suitable for ggplot/ggtree
##
## Suggested Enhancements:
##   - Batch extraction for multiple time points
##   - Add branch width scaling by rate or HPD size
##   - Map projection support
##
## @param timePt
##   Numeric: decimal year (or time unit) at which to slice the tree
##
## @param tr
##   MCC tree object containing:
##       - tr$edge, tr$latlon, tr$nodeTimes
##       - tr$hpds / tr$hpd_diffs for HPD polygons
##       - Optional: tr$fullset for discrete trait interpolation
##
## @param ii
##   Optional numeric vector of branch indices. Auto-selected if empty.
##
## @param xlim
##   Numeric vector length 2: longitude limits
##
## @param ylim
##   Numeric vector length 2: latitude limits
##
## @param excludedTips
##   Character vector of tip labels to exclude from interpolation
##
## @param colorBy
##   "none" (default) | "HPDwidth" | "trait": used for coloring HPDs
##
## @param minHPDwidth
##   Minimum HPD width threshold (only branches above included)
##
## @return
##   List containing:
##       - pts: interpolated x/y coordinates
##       - pts_polys: HPD polygons
##       - ii: branch indices used
##       - fractTime: fraction along branch where timePt falls
##       - Optional:
##           - interpol_trait: interpolated discrete traits
##           - midSet: interpolated trait matrices
##       - df: ggplot-ready data frame with x, y, branch, HPD, trait info
############################################################

pts_and_hpds_at_time <- function(timePt,
                                 tr,
                                 ii=c(),
                                 xlim=c(-180,180),
                                 ylim=c(-60,80),
                                 excludedTips=c(),
                                 colorBy="none",
                                 minHPDwidth=0) {

  stopifnot(!is.null(tr$edge), !is.null(tr$latlon), !is.null(tr$nodeTimes))
  
  # Extract from/to coordinates
  fromY <- tr$latlon[tr$edge[,1],1]
  fromX <- tr$latlon[tr$edge[,1],2]
  toY   <- tr$latlon[tr$edge[,2],1]
  toX   <- tr$latlon[tr$edge[,2],2]
  
  fromTime <- tr$nodeTimes[tr$edge[,1]]
  toTime   <- tr$nodeTimes[tr$edge[,2]]
  
  branchTime <- toTime - fromTime
  branchX    <- toX - fromX
  branchY    <- toY - fromY
  
  # HPDs
  hpds      <- tr$hpds
  hpd_diffs <- tr$hpd_diffs
  
  # Auto-select branches intersecting timePt and within limits
  if(length(ii)==0){
    ii <- which(
      toTime >= timePt & fromTime <= timePt &
      toX >= xlim[1] & toX <= xlim[2] &
      toY >= ylim[1] & toY <= ylim[2]
    )
  }
  
  # Exclude specified tips
  if(length(excludedTips) > 0){
    ex_tips <- match(excludedTips, tr$tip.label)
    ex_tips <- ex_tips[is.finite(ex_tips)]
    if(length(ex_tips) > 0){
      ee <- match(ex_tips, tr$edge[,2])
      ii <- setdiff(ii, ee)
    }
  }
  
  # Compute fraction along branch
  fractTime <- 1 - (toTime - timePt) / branchTime
  
  # Filter branches by minimum HPD width if requested
  if(minHPDwidth > 0 && length(ii) > 0){
    widths <- sapply(ii, function(j) {
      max(hpd_diffs[[j]]) - min(hpd_diffs[[j]])
    })
    ii <- ii[widths >= minHPDwidth]
  }
  
  # Initialize outputs
  if(length(ii)==0) return(list(pts=NULL, pts_polys=NULL, ii=NULL, fractTime=NULL))
  if(length(ii)==1) ii <- c(ii, ii)
  
  # Interpolate coordinates
  x_pos <- (fractTime * branchX + fromX)[ii]
  y_pos <- (fractTime * branchY + fromY)[ii]
  pts <- list(x=x_pos, y=y_pos)
  
  # Interpolate HPDs
  pts_polys <- lapply(seq_along(ii), function(j){
    fromNode <- tr$edge[ii[j],1]
    fractTime[ii[j]] * hpd_diffs[[ii[j]]] + hpds[[fromNode]]
  })
  
  # Interpolate discrete traits if available
  interpol_trait <- midSet <- NULL
  if(!is.null(tr$fullset)){
    fromSet <- tr$fullset[tr$edge[ii,1],]
    toSet   <- tr$fullset[tr$edge[ii,2],]
    midSet  <- interpolate_discrete_element(fromSet, toSet, fractTime[ii])
    interpol_trait <- colnames(midSet)[apply(midSet,1,which.max)]
  }
  
  # Prepare ggplot-ready data frame
  df <- data.frame(
    branch=ii,
    x=x_pos,
    y=y_pos,
    fractTime=fractTime[ii]
  )
  if(!is.null(interpol_trait)) df$trait <- interpol_trait
  
  # Color coding (optional)
  if(colorBy=="HPDwidth"){
    df$HPDwidth <- sapply(ii, function(j) max(hpd_diffs[[j]]) - min(hpd_diffs[[j]]))
  } else if(colorBy=="trait" && !is.null(interpol_trait)){
    df$HPDwidth <- NULL
    df$traitCol <- interpol_trait
  }
  
  return(list(
    pts=pts,
    pts_polys=pts_polys,
    ii=ii,
    fractTime=fractTime[ii],
    interpol_trait=interpol_trait,
    midSet=midSet,
    df=df
  ))
}
