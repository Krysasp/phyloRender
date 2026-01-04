get_BEAST_cols <- function(n, sat=0.6, bright=0.8, transparency=1) {
  grDevices::hsv(
    h = seq(0, 1, length.out = n + 1)[-1],
    s = sat,
    v = bright,
    alpha = transparency
  )
}

colourize_BEAST_cols <- function(states) {
  ustates <- sort(unique(states))
  ucols   <- get_BEAST_cols(length(ustates))
  statecols <- ucols[match(states, ustates)]
  list(ustates=ustates, ucols=ucols, statecols=statecols)
}