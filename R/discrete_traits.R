############################################################
## SCRIPT CREDENTIALS 
############################################################
## Tool       :  # phyloRender as part of the 'phyloAct' program
## Current script    :  # discrete_traits.R
## Author            :  # Jonathan Chan
## Affiliation       :  # Institute of Health and Community Medicine
## Last Updated      :  # 2025-12-27
############################################################
## FUNCTION: addDiscreteTraits / addFullTraitSet / addContinuousTraits
############################################################
## Description:
##   Extracts discrete trait probabilities from BEAST MCC annotations
##   and stores them in a matrix form for plotting and interpolation.
##   Optionally, continuous trait values (e.g., substitution rates, branch lengths)
##   can be extracted for each node and tip using addContinuousTraits().
##
## Features:
##   - Compatible with main.R plotting functions
##   - Prepares discrete trait matrix for tip & internal node interpolation
##   - Allows adding continuous traits to tree object
##   - Useful for animating discrete/continuous trait evolution over time
##
## Suggested Enhancements:
##   - Color-code traits automatically for plotting
##   - Use with ggtree or plotly for interactive visualization
##   - Scale branch width by continuous trait values
##
## @param tr
##   MCC tree object with tr$details annotations
##
## @param propIndex
##   Integer specifying which trait column to extract / use
##
## @param continuousKeys
##   Character vector of annotation keys to extract as continuous traits
##   Example: c("rate","length","location.rate")
##
## @return
##   Updated tree object with:
##     - tr$props       : discrete trait matrix
##     - tr$uprops      : unique values per discrete trait
##     - tr$propNames   : discrete trait names
##     - tr$fullset     : full discrete probability matrix for internal nodes
##     - tr$contTraits  : named matrix of continuous trait values per node/tip
############################################################

addDiscreteTraits <- function(tr) {

  # Regex pattern to detect discrete traits from BEAST annotation strings
  discreteRegex <- "[\\,\\&][A-Za-z0-9_\\-]+\\.prob"
  pos <- gregexpr(discreteRegex, tr$details[1])[[1]]

  if (pos[1] < 0) stop("No discrete traits found in tree annotations.")

  # Extract trait names
  is <- pos + 1
  ie <- is + attributes(pos)$match.length - 7
  traitNames <- substring(tr$details[1], is, ie)
  numTraits <- length(traitNames)

  # Initialize trait matrices
  props <- matrix(NA, nrow=length(tr$details), ncol=numTraits)
  colnames(props) <- traitNames
  uprops <- vector("list", numTraits)

  # Fill in trait probabilities for each node/tip
  for (j in seq_len(numTraits)) {
    vals <- apply(as.matrix(tr$details), 1, function(x) {
      v <- sub(paste0(".*",traitNames[j],"="),"",x)
      v <- strsplit(v,",")[[1]][1]
      gsub("[\\]\\\"]","",v)
    })
    props[,j] <- vals
    uprops[[j]] <- sort(unique(vals))
  }

  tr$props     <- props
  tr$uprops    <- uprops
  tr$propNames <- traitNames

  return(tr)
}

# Full discrete trait matrix for internal nodes & tips
getFullTraitSet <- function(tr, propIndex=1) {
  utraits <- tr$uprops[[propIndex]]
  trait   <- tr$propNames[propIndex]

  ntips <- length(tr$tip.label)
  ancs  <- (ntips+1):(ntips+tr$Nnode)

  res <- lapply(tr$details[ancs], getTraitSet,
                traitName=trait,
                orderTraits=FALSE,
                utraits=utraits)

  mat <- do.call(rbind, res)

  tipmat <- matrix(0, nrow=ntips, ncol=length(utraits))
  colnames(tipmat) <- utraits
  for (i in seq_len(ntips)) {
    tipmat[i, match(tr$props[i,propIndex], utraits)] <- 1
  }

  fullset <- rbind(tipmat, mat)
  colnames(fullset) <- utraits
  return(fullset)
}

# Add the full probability set to tree object
addFullTraitSet <- function(tr, propIndex=1) {
  tr$fullset <- getFullTraitSet(tr, propIndex)
  return(tr)
}

# Interpolate between discrete traits along branches
interpolate_discrete_element <- function(fromSet, toSet, fractTime) {
  ff <- matrix(fractTime, nrow=length(fractTime), ncol=ncol(fromSet))
  (fromSet * (1-ff)) + (toSet * ff)
}

############################################################
## New Function: addContinuousTraits
############################################################
## Description:
##   Extracts continuous trait values from BEAST MCC annotations
##   and stores them in tr$contTraits.
##   Works with values such as rate, branch length, location.rate, etc.
##
## @param tr
##   MCC tree object with tr$details annotations
##
## @param continuousKeys
##   Character vector of annotation keys to extract, e.g., c("rate","length")
##
## @return
##   Updated tree object with tr$contTraits: numeric matrix [nodes x traits]
############################################################

addContinuousTraits <- function(tr, continuousKeys=c("rate","length")) {

  if(is.null(tr$details)) stop("Tree object lacks $details annotations.")

  nNodes <- length(tr$details)
  contMat <- matrix(NA, nrow=nNodes, ncol=length(continuousKeys))
  colnames(contMat) <- continuousKeys

  # Loop through keys and parse numeric values
  for(k in seq_along(continuousKeys)) {
    key <- continuousKeys[k]
    contMat[,k] <- sapply(tr$details, function(x) {
      m <- regmatches(x, regexpr(paste0(key,"=\\{?[0-9Ee\\.\\-]+"), x))
      if(length(m)==0) return(NA)
      as.numeric(sub(paste0(key,"="),"",m))
    })
  }

  tr$contTraits <- contMat
  return(tr)
}