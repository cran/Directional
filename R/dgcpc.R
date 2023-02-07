dgcpc <- function(x, omega, g, rho, rads = FALSE, logden = FALSE) {
  if ( !rads )  x <- x * pi/180
  x <- cbind( cos(x), sin(x) )
  ksi <- c( cos(omega), sin(omega) )
  mu <- g * ksi
  sinv <- matrix( c( ksi[1]^2 + ksi[2]^2/rho, ksi[1] * ksi[2] * (1 - 1/rho),
                  ksi[1] * ksi[2] * (1 - 1/rho), ksi[2]^2 + ksi[1]^2/rho ), ncol = 2 )
  a <- as.vector( x %*% mu )
  b <- Rfast::rowsums( x %*% sinv * x)
  den <-  - 0.5 * log(rho) - log(2 * pi) - log( b * sqrt(g^2 + 1) - a * sqrt(b) )
  if ( !logden )  den <- exp(den)
  den
}
