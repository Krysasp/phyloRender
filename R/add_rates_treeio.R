############################################################
## SCRIPT CREDENTIALS 
############################################################
## Tool       :  # phyloRender as part of the 'phyloAct' program
## Current script    :  # read_beast_mcc.R
## Author            :  # Jonathan Chan
## Affiliation       :  # Institute of Health and Community Medicine
## Last Updated      :  # 2025-12-05
############################################################
## FUNCTION: addRates_treeio
############################################################
## Description:
##   Extracts substitution rate estimates and 95% HPD intervals
##   from treeio MCC tree annotations and adds them to the tree object.
##   Works with trees that have a metadata column `details` containing
##   rate annotations from BEAST runs.
##
## Parameters:
##   tr         : A treeio object (e.g., phylo or treedata object) that
##                must contain a 'details' column with branch rate
##                annotations from BEAST.
##
## Usage Example:
##   tr <- addRates_treeio(tr)
##
## Parameter Modifications / Customization:
##   - patt: You can modify the pattern used to extract HPD columns,
##           e.g., change "rate_95%HPD_1" or "rate_95%HPD_2" if
##           your BEAST output uses different labels.
##
## Suggested Enhancements / New Functionalities:
##   1. Automatically calculate median rates from HPD intervals.
##   2. Add coloring metadata for branches by rate for visualization.
##   3. Add option to compute log-scaled rates if rate distribution is skewed.
##   4. Integrate with animation-ready trees for gganimate (time slices).
##   5. Provide a warning if HPD extraction fails for any branch.
##
## Notes:
##   - Requires treeio object with a 'details' column
##   - Currently extracts numeric values from BEAST-style annotation strings
############################################################

addRates_treeio <- function(tr) {
  
  # 1. Ensure tree has 'details' metadata
  stopifnot(!is.null(tr$details))
  
  # 2. Helper function to extract numeric values by pattern
  extract <- function(x, patt) {
    # Extract the first matching pattern in the comma-separated string
    as.numeric(sub(paste0(".*", patt, "="), "", strsplit(x, ",")[[1]][1]))
  }
  
  # 3. Extract lower and upper 95% HPD for substitution rates
  tr$rate_range_1 <- sapply(tr$details, extract, patt = "rate_95%HPD_1")
  tr$rate_range_2 <- sapply(tr$details, extract, patt = "rate_95%HPD_2")
  
  # 4. Optional: calculate median rate
  tr$rate_median <- (tr$rate_range_1 + tr$rate_range_2) / 2
  
  # 5. Optional: add log-scaled rates for skewed distributions
  tr$rate_log <- log(tr$rate_median + 1e-8) # avoid log(0)
  
  # 6. Optional: add rate-based color categories for plotting
  tr$rate_category <- cut(tr$rate_median,
                          breaks = quantile(tr$rate_median, probs = seq(0, 1, 0.25), na.rm = TRUE),
                          include.lowest = TRUE,
                          labels = c("low", "medium-low", "medium-high", "high"))
  
  return(tr)
}

############################################################
## Required supporting packages/functions:
## - treeio: for reading MCC trees
## - ape: optional if converting tree to phylo for plotting
## - ggtree: for visualization and color-mapped rates
############################################################
