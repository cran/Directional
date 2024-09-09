# Density function
# @param y, vector with length d
# @param B, vector with length (d+2)(d-1)/2
# @return density of ESAG with parameter B = c(mu,gamma) evaluated at y.
dESAGd <- function(y, mu, gam, logden = FALSE) {

  y <- as.matrix(y)
  if  ( dim(y)[2] == 1 )   y <- t(y)
  n <- dim(y)[1]  ;  d <- dim(y)[2]
  R <- .rotation(gam)
  lambda <- .parameter(gam)[[ 1 ]]

  eigenvector_hat <- .ONB(mu)
  eigenvector <- eigenvector_hat[, -d] %*% R
  P <- cbind(eigenvector, eigenvector_hat[, d] )
  V <- P %*% ( t(P) * c(lambda, 1) )

  p <- d - 1
  yVy <- Rfast::mahala(y, numeric(d), V)
  Cd <- (2 * pi)^(-0.5 * p)
  a <- as.vector( y %*% mu / sqrt(yVy) )
  Mp <- matrix(1, nrow = n, ncol = p)
  Mp[, 1] <- a * pnorm(a) + dnorm(a)
  Mp[, 2] <- (1 + a^2) * pnorm(a) + a * dnorm(a)
  if ( p >= 3 )   for ( j in 3:p )  Mp[, j] <- a * Mp[, j - 1] + (j - 1) * Mp[, j - 2]
  l <- log(Cd) - 0.5 * d * log(yVy) + 0.5*( a^2 - sum(mu^2) ) + log( Mp[, p] )
  if ( logden )  {
    return( l )
  } else   return( exp(l) )
}



# dESAG_old <- function(y, mu, gamma) {
#   d <- length(y)
#   R <- .rotation(gamma)
#   lambda <- .parameter(gamma)[[ 1 ]]
#   V_inv <- .covariance_inv_matrix(mu, lambda, R)
#   p <- d - 1
#   yVy <- t(y) %*% V_inv %*% y
#   Cd <- (2 * pi)^(-0.5 * p)
#   a <- as.vector( t(y) %*% mu / sqrt(yVy) )
#   Mp <- numeric(p)
#   Mp[1] <- a * pnorm(a) + dnorm(a)
#   Mp[2] <- (1 + a^2) * pnorm(a) + a * dnorm(a)
#   if ( p >= 3 )   for ( i in 3:p )  Mp[i] <- a * Mp[i - 1] + (i - 1) * Mp[i - 2]
#   M_p <- Mp[p]
#   Cd/( yVy^(d/2) ) * exp(0.5*( ( t(y) %*% mu)^2 / yVy - t(mu) %*% mu )) * M_p
# }
