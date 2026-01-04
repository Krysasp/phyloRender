getEl <- function(x, sep, ind=1, fromEnd=FALSE,
                  ex=0, reconstruct=FALSE) {

  if (is.na(x)) return(NA)

  sp <- strsplit(x, sep)[[1]]

  if (fromEnd) {
    ind <- length(sp) - ind + 1
  }

  if (reconstruct) {
    if (ex > 0) {
      return(paste(sp[-(1:ex)], collapse=sep))
    }
  }

  return(sp[ind])
}