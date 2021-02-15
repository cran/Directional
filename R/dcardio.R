dcardio <- function(x, m, rho, rads = FALSE, logden = FALSE) {
  if ( !rads )  x <- x * pi/180
  den <-  -log(2 * pi) + log1p( 2 * rho * cos(x - m) )
  if ( !logden )  den <- exp(den)
  den
}
