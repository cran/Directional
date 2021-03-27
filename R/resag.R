resag <- function(n, mu, gam) {
  ## n is the sample size
  ## para contains the mean vector and the two gammas
  gam1 <- gam[1]
  gam2 <- gam[2]
  if ( gam1 == 0  &  gam2 == 0 ) {
    y <- Rfast::matrnorm(n, 3)  ## IAG is used
    y <- Rfast::eachrow(y, mu, oper = "+")
  } else {
    m0 <- sqrt( mu[2]^2 + mu[3]^2 )
    rl <- sqrt( sum( mu^2 ) )
    x1b <- c( -m0^2, mu[1] * mu[2], mu[1] * mu[3] ) / (m0 * rl)
    x2b <- c( 0, -mu[3], mu[2] )/m0
    T1 <- tcrossprod( x1b )
    T2 <- tcrossprod( x2b )
    T12 <- tcrossprod( x1b, x2b )
    vinv <- diag(3) + gam1 * ( T1 - T2 ) + gam2 * ( T12 + t(T12) ) + ( sqrt(gam1^2 + gam2^2 + 1) - 1 ) * ( T1 + T2 )
    va <- solve(vinv)
    y <- Rfast::rmvnorm(n, mu, va)  ## sample from a multivariate normal in R^3
  }
  colnames(y) <- c("X", "Y", "Z")
  y / sqrt( Rfast::rowsums(y^2) )  ## ESAG simulated values
}
