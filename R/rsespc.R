rsespc <- function(n, mu, theta) {
  the1 <- theta[1]
  the2 <- theta[2]
  if ( the1 == 0  &  the2 == 0 ) {
    y <- Rfast::rmvt(n, mu, sigma = diag(3), v = 1)
  } else {
    m0 <- sqrt(mu[2]^2 + mu[3]^2)
    rl <- sqrt(sum(mu^2))
    x1b <- c(-m0^2, mu[1] * mu[2], mu[1] * mu[3])/(m0 * rl)
    x2b <- c(0, -mu[3], mu[2])/m0
    T1 <- tcrossprod(x1b)
    T2 <- tcrossprod(x2b)
    T12 <- tcrossprod(x1b, x2b)
    vinv <- diag(3) + the1 * (T1 - T2) + the2 * (T12 + t(T12)) + 
            ( sqrt(the1^2 + the2^2 + 1) - 1 ) * (T1 + T2)
    va <- solve(vinv)
    y <- Rfast::rmvt(n, mu, va, v = 1)
  }  
  colnames(y) <- c("X", "Y", "Z")
  y / sqrt( Rfast::rowsums(y^2) )
}
