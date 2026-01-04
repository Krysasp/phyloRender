############################################################
## SCRIPT CREDENTIALS 
############################################################
## Tool       :  # phyloRender as part of the 'phyloAct' program
## Current script    :  # read_beast_mcc.R
## Author            :  # Jonathan Chan
## Affiliation       :  # Institute of Health and Community Medicine
## Last Updated      :  # 2025-12-05
############################################################
## FUNCTION: read_mcc_treeio
############################################################
## Description:
##   Reads a BEAST MCC tree using the `treeio` package with a 
##   fallback to a text-based parser if treeio fails.  
##   Extracts branch annotations and node heights into a 
##   phylo object with additional metadata for downstream analysis.
##
## Parameters:
##   file           : Character string. Path to the BEAST MCC tree file (NEXUS or Newick format).
##   prefer_treeio  : Logical (TRUE/FALSE). If TRUE, attempt to read the tree using treeio first; 
##                    if FALSE, skip treeio and go directly to text-based parsing.
##
## Return Value:
##   A phylo object with additional metadata:
##     - tr$details    : BEAST annotations per node (rates, traits, etc.)
##     - tr$nodeTimes  : Node ages relative to tree root (if height available)
##
## Usage Example:
##   tr <- read_mcc_treeio("example_mcc.nexus", prefer_treeio=TRUE)
##
## Parameter Modifications / Customization:
##   - file: change this to the path of your BEAST MCC output file.
##   - prefer_treeio: set FALSE if treeio fails or if you prefer a pure text parser.
##   - If your BEAST output has non-standard annotation names, you may need to adjust 
##     how `tr$details` or `tr$nodeTimes` are calculated.
##
## Suggested Enhancements / New Functionalities:
##   1. Automatically detect BEAST1 vs BEAST2 output and adjust parsing.
##   2. Extract discrete traits or spatial metadata automatically.
##   3. Add warning if treeio parsing succeeds but some nodes lack annotations.
##   4. Integrate branch rate and HPD extraction for immediate visualization.
##   5. Optionally return a ggtree-ready object with mapped annotations for plotting.
##
## Notes:
##   - Requires `treeio` package for treeio-based parsing.
##   - Requires fallback function `read_mcc_tr()` for text-based parsing.
##   - Works best with trees exported from BEAST as MCC (maximum clade credibility) trees.
############################################################

read_mcc_treeio <- function(file, prefer_treeio=TRUE) {

  if (prefer_treeio) {
    # Attempt treeio-based reading
    res <- try(treeio::read.beast(file), silent=TRUE)
    if (!inherits(res, "try-error")) {

      # Convert to phylo object
      tr <- treeio::as.phylo(res)
      tr$node.label <- NULL

      # Extract annotations and order by node
      ann <- res@data
      ann <- ann[order(ann$node),]

      tr$details <- ann$annotation

      # Calculate node times relative to tree root if available
      if ("height" %in% names(ann)) {
        tr$nodeTimes <- max(ann$height) - ann$height
      }

      return(tr)
    }
  }

  # Fallback to text-based MCC reader
  message("Falling back to text-based MCC reader")
  read_mcc_tr(file, BEAST2=TRUE)
}
