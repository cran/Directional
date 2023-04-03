rgcpc <- function(n, mu = NULL, omega, g, rho, rads = TRUE) {
  if ( is.null(mu) ) {
    if ( !rads )  omega <- omega * pi/180
    ksi <- cbind( cos(omega), sin(omega) )
    mu <- g * ksi
  } else  ksi <- mu / sqrt( sum(mu^2) )
  sinv <- matrix( c( ksi[1]^2 + ksi[2]^2/rho, ksi[1] * ksi[2] * (1 - 1/rho), 
                  ksi[1] * ksi[2] * (1 - 1/rho), ksi[2]^2 + ksi[1]^2/rho ), ncol = 2 )
  s <- solve(sinv)
  x <- Rfast::rmvt(n, mu, s, 1)
  x <- x / sqrt( Rfast::rowsums(x^2) )
  ( atan(x[, 2]/x[, 1]) + pi * I(x[, 1] < 0) ) %% (2 * pi)
}