dcircexp <- function(x, lambda, rads = FALSE, logden = FALSE) {
  if ( !rads )  x <- x * pi/180
  den <- log(lambda) - lambda * x - log( 1 - exp(-2 * pi * lambda) )
  if ( !logden )  den <- exp(den)
  den
}
