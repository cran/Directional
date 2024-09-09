iagd <- function(y, mu, logden = FALSE) {
  y <- as.matrix(y)
  d <- dim(y)[2]
  if ( d == 1 )   y <- t(y)
  n <- dim(y)[1]  ;  d <- dim(y)[2]

  if ( d == 3 ) {
    a <- tcrossprod(mu, y)
    a2 <- a^2
    pa <- pnorm(a)
    da <- dnorm(a)
    gm <- pa + a2 * pa + a * da
    den <-  -log(2 * pi) + 0.5 * a2 - 0.5 * sum(mu^2) + log(gm)

  } else if ( d > 3 ) {
    p <- d - 1
    Cd <- (2 * pi)^(-0.5 * p)
    a <- as.vector( y %*% mu )
    Mp <- matrix(1, nrow = n, ncol = p)
    Mp[, 1] <- a * pnorm(a) + dnorm(a)
    Mp[, 2] <- (1 + a^2) * pnorm(a) + a * dnorm(a)
    if ( p >= 3 )   for ( j in 3:p )  Mp[, j] <- a * Mp[, j - 1] + (j - 1) * Mp[, j - 2]
    den <- log(Cd) + 0.5*( a^2 - sum(mu^2) ) + log( Mp[, p] )
  }

  if ( !logden )  den <- exp(den)
  den
}



.IAGd <- function(y, mu, logden = FALSE) {
  y <- as.matrix(y)
  if  ( dim(y)[2] == 1 )   y <- t(y)
  n <- dim(y)[1]  ;  d <- dim(y)[2]
  p <- d - 1
  Cd <- (2 * pi)^(-0.5 * p)
  a <- as.vector( y %*% mu )
  Mp <- matrix(1, nrow = n, ncol = p)
  Mp[, 1] <- a * pnorm(a) + dnorm(a)
  Mp[, 2] <- (1 + a^2) * pnorm(a) + a * dnorm(a)
  if ( p >= 3 )   for ( j in 3:p )  Mp[, j] <- a * Mp[, j - 1] + (j - 1) * Mp[, j - 2]
  l <- log(Cd) + 0.5*( a^2 - sum(mu^2) ) + log( Mp[, p] )
  if ( logden )  {
    return( l )
  } else   return( exp(l) )
}
