dwrapcauchy <- function(x, m, rho, rads = FALSE, logden = FALSE) {
  if ( !rads )  x <- x * pi/180
  den <-  -log(2 * pi) + log(1 - rho^2) - log1p( rho^2 - 2 * rho * cos(x - m) )
  if ( !logden )  den <- exp(den)
  den
}

