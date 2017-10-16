ESAGdensity <- function(y, param, logden = FALSE) {
  ## y is a matrix withe the sphericla data, unit vectors
  ## m is the mean vector
  ## param contans a) the mean vector and the two gammas
  gam1 <- param[4]
  gam2 <- param[5]
  m <- param[1:3]
  m0 <- sqrt( m[2]^2 + m[3]^2 )
  rl <- sqrt( sum( m^2 ) )
  x1b <- c( -m0^2, m[1] * m[2], m[1] * m[3] ) / (m0 * rl)
  x2b <- c( 0, -m[3], m[2] )/m0
  T1 <- tcrossprod( x1b )
  T2 <- tcrossprod( x2b )
  T12 <- tcrossprod( x1b, x2b )
  va <- diag(3) + gam1 * ( T1 - T2 ) + gam2 * ( T12 + t(T12) ) + ( sqrt(gam1^2 + gam2^2 + 1) - 1 ) * ( T1 + T2 )
  y <- as.matrix(y)
  if   ( dim(y)[2] == 1 )   y <- t(y)
  g2 <- as.vector( y %*% m )
  g1 <- Rfast::rowsums( y %*% va * y )
  a <- g2 / sqrt(g1)
  M2 <- ( 1 + a^2 ) * pnorm(a) + a * dnorm(a)
  l <-  - log(2 * pi) + 0.5 * a^2 - 0.5 * sum(m^2) - 1.5 * log(g1) + log(M2)
  if ( logden )  {
    return( l )
  } else   return( exp(l) )
}
