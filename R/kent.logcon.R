## Normalising constant of the Kent distribution

kent.logcon <- function(k, b, j = 100) {
  j <- 0:j
  ka <- 2 * pi * gamma(j + 0.5) / gamma(j + 1)* b^( 2 * j ) * ( k / 2 )^( -2 * j - 0.5 ) * besselI(k, 2 * j + 0.5)
  log( sum(ka) )
}
