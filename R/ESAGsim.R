ESAGsim <- function(n, param) {
  ## n is the sample size
  ## para contains the mean vector and the two gammas
  m <- param[1:3]
  gam1 <- param[4]
  gam2 <- param[5]
  if ( gam1 == 0  &  gam2 == 0 ) {
    y <- Rfast::matrnorm(n, 3)  ## IAG is used
    y <- Rfast::eachrow(y, m, oper = "+")
  } else {
    m0 <- sqrt( m[2]^2 + m[3]^2 )
    rl <- sqrt( sum( m^2 ) )
    x1b <- c( -m0^2, m[1] * m[2], m[1] * m[3] ) / (m0 * rl)
    x2b <- c( 0, -m[3], m[2] )/m0
    T1 <- tcrossprod( x1b )
    T2 <- tcrossprod( x2b )
    T12 <- tcrossprod( x1b, x2b )
    vinv <- diag(3) + gam1 * ( T1 - T2 ) + gam2 * ( T12 + t(T12) ) + ( sqrt(gam1^2 + gam2^2 + 1) - 1 ) * ( T1 + T2 )
    va <- solve(vinv)
    y <- Rfast::rmvnorm(n, m, va)  ## sample from a multivariate normal in R^3
  }
  colnames(y) <- c("X", "Y", "Z")
  y / sqrt( Rfast::rowsums(y^2) )  ## ESAG simulated values
}
