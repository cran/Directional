rcipc <- function(n, mu = NULL, omega, g, rads = TRUE) {
  if ( is.null(mu) ) {
    if ( !rads )  omega <- omega * pi/180
    mu <- g * c( cos(omega), sin(omega) )
  } 
  x <- Rfast::rmvt(n, mu, diag(2), 1)
  x <- x / sqrt( Rfast::rowsums(x^2) )
  ( atan(x[, 2]/x[, 1]) + pi * I(x[, 1] < 0) ) %% (2 * pi)
}  