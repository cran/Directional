## Probability density function of the
## von Mises-Fisher distribution
## May 2016
## References: Arthur Pewsey, Markus Neuhauser, and Graeme D. Ruxton (2013)
## Circular Statistics in R
pvm <- function(u, m, k, rads = FALSE) {
  if ( !rads )  {
     u <- u * pi / 180
     m <- m * pi / 180
  }
  if ( k > 0 ) {
    f <- 2 * pi * besselI(k, 0)
    funa <- function(u, m, k)  exp(k * cos(u - m) )
    n <- length(u)
    lo <- rep(0, n)
    prob <- mapply( function(lo, u, m, k)  integrate(funa, lo, u, m, k)$value, lo, u, m = m, k = k) / f
  } else  prob <- u / ( 2 * pi )
  prob
}
