pcircexp <- function(u, lambda, rads = FALSE) {
  if ( !rads )  u <- u * pi / 180
  ( 1 - exp(-lambda * u) ) / ( 1 - exp(-2 * pi * lambda) )
}
