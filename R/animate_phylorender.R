############################################################
## Tool       :  # phyloRender as part of the 'phyloAct' program
## Current script    :  # animate_phyloRender.R
## Author            :  # Jonathan Chan
## Affiliation       :  # Institute of Health and Community Medicine
## Last Updated      :  # 2025-10-28
############################################################
## FUNCTION: animate_tree_spatial / animate_phyloRender
############################################################
## Description:
##   Creates a series of frames showing the spatial phylodynamics 
##   of a MCC tree through time.
##   Automatically extracts interpolated tip positions and HPDs
##   at each time point using pts_and_hpds_at_time().
##   Can optionally show discrete trait interpolation and color-coding.
##
## Features:
##   - Saves PNG frames for each time point in `outdir`
##   - Plots HPD polygons (semi-transparent) along branches
##   - Plots tip positions with optional color/trait annotation
##   - Supports variable spatial limits (`xlim`, `ylim`)
##   - Optional trait coloring if tr$fullset exists
##   - Compatible with main.R plotting workflow
##
## Suggested Enhancements:
##   - Direct gif/mp4 export using magick/gganimate
##   - Add branch width scaling based on rate or HPD size
##   - Support interactive plots using ggtree or plotly
##
## @param tr
##   MCC tree object containing:
##       - tr$edge, tr$latlon, tr$nodeTimes
##       - tr$hpds / tr$hpd_diffs for HPD polygons
##       - Optional: tr$fullset for discrete trait interpolation
##
## @param times
##   Numeric vector: time points (decimal years) at which to slice tree
##
## @param outdir
##   Character: output directory where frames will be saved
##
## @param prefix
##   Character: frame file prefix
##
## @param xlim
##   Numeric vector length 2: longitude limits of plotting area
##
## @param ylim
##   Numeric vector length 2: latitude limits of plotting area
##
## @param propIndex
##   Integer: which discrete trait column to show, 0 = none
##
## @param colorBy
##   Character: "none" | "HPDwidth" | "trait" (colors tip/HPD accordingly)
##
## @param minHPDwidth
##   Numeric: only branches with HPD width above this threshold are included
############################################################

animate_phyloRender <- function(
  tr,
  times,
  outdir="frames",
  prefix="frame",
  xlim=c(-180,180),
  ylim=c(-60,80),
  propIndex=0,
  colorBy="none",
  minHPDwidth=0
) {

  # Create output folder if not exists
  if (!dir.exists(outdir))
    dir.create(outdir, recursive=TRUE)
  
  # Loop through each time point
  for(i in seq_along(times)) {
    
    # Prepare PNG device
    png(
      file.path(outdir, sprintf("%s_%04d.png", prefix, i)),
      width=1200, height=800
    )
    
    # Plot MCC tree as background
    plot_mcc_tree_with_hpds(tr, xlim=xlim, ylim=ylim, propIndex=propIndex)
    
    # Extract interpolated positions and HPDs at current time
    slice <- pts_and_hpds_at_time(
      timePt = times[i],
      tr = tr,
      xlim = xlim,
      ylim = ylim,
      colorBy = colorBy,
      minHPDwidth = minHPDwidth
    )
    
    # Plot HPD polygons
    if(!is.null(slice$pts_polys)) {
      for(j in seq_along(slice$pts_polys)) {
        colHPD <- if(colorBy=="HPDwidth" && !is.null(slice$df$HPDwidth)){
          rgb(1,0,0, alpha=min(0.8, slice$df$HPDwidth[j]/max(slice$df$HPDwidth)))
        } else if(colorBy=="trait" && !is.null(slice$df$traitCol)){
          grDevices::adjustcolor(slice$df$traitCol[j], alpha.f=0.4)
        } else {
          rgb(1,0,0,0.3)
        }
        graphics::polygon(slice$pts_polys[[j]], border=NA, col=colHPD)
      }
    }
    
    # Plot tip points
    if(!is.null(slice$pts)) {
      ptCols <- if(colorBy=="trait" && !is.null(slice$df$traitCol)){
        slice$df$traitCol
      } else {
        "red"
      }
      graphics::points(slice$pts$x, slice$pts$y, pch=21, bg=ptCols, cex=1.2)
    }
    
    # Add title showing time
    graphics::title(paste("Time:", round(times[i],3)))
    
    # Close PNG device
    dev.off()
  }
}
