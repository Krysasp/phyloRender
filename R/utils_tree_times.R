distFromRoot <- function(tr, node) {
  ntips <- length(tr$tip.label)
  root  <- ntips + 1

  if (node == root) return(0)

  parent_edge <- which(tr$edge[,2] == node)
  parent      <- tr$edge[parent_edge,1]
  bl          <- tr$edge.length[parent_edge]

  return(bl + distFromRoot(tr, parent))
}

nodeTimes <- function(tr, youngestTip) {

  ntips <- length(tr$tip.label)
  nnodes <- ntips + tr$Nnode
  nodeTimes <- numeric(nnodes)

  for (i in 1:nnodes) {
    nodeTimes[i] <- youngestTip - distFromRoot(tr, i)
  }

  tr$nodeTimes <- nodeTimes
  return(tr)
}