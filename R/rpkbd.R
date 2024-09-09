rpkbd <- function(n, mu, rho) {

  lam <- 2 * rho / (1 + rho^2)
  d <- length(mu)

  cdlb <- function(be, lam, d) {
    d * log( ( 1 + sqrt(1 - lam^2) ) / (1 + sqrt(1 - lam^2 / be) ) ) - log( 1 - be)
  }
  bstar <- optimize(cdlb, c(lam * (2 - lam), 1), lam = lam, d = d)$minimum
  b1 <- bstar / (1 - bstar)
  b2 <-  - 1 + 1 / sqrt(1 - bstar)

  muv <- as.vector(mu)
  u <- Rfast2::Runif(2 * n)
  z <- Rfast::matrnorm(2 * n, d)
  mz <- as.vector(z %*% muv)
  ##zz <- Rfast::mahala(z, numeric(d), diag(d) )
  zz <- Rfast::rowsums(z^2)
  com <- sqrt( zz + b1 * mz^2 )
  qa <- ( mz * (1 + b2) ) / com
  wa <- log(u) <= 0.5 * d * ( - log(1 - lam * qa) + log(1 - bstar * qa^2) -
                                log(2 / (1 + sqrt(1 - lam^2/bstar) ) ) )
  wa <- which(wa)
  n1 <- length(wa)
  if ( n1 > 0 ) {
    mu <- matrix( rep(muv, each = n1), ncol = d)
    x <- (z[wa, ] + b2 * mz[wa] * mu) / com[wa]
  } else x <- NULL

  while (n1 < n) {
    u <- Rfast2::Runif(2 *(n - n1) )
    z <- Rfast::matrnorm(2 * (n - n1 ), d)
    mz <- as.vector(z %*% muv)
    ##zz <- Rfast::mahala(z, numeric(d), diag(d) )
	zz <- Rfast::rowsums(z^2)
    com <- sqrt( zz + b1 * mz^2 )
    qa <- ( mz * (1 + b2) ) / com
    wa <- log(u) <= 0.5 * d * ( - log(1 - lam * qa) + log(1 - bstar * qa^2) -
                                  log(2 / (1 + sqrt(1 - lam^2/bstar) ) ) )
    wa <- which(wa)
    n2 <- length(wa)
    if ( n2 > 0 ) {
      mu <- matrix( rep(muv, each = n2), ncol = d)
      x <- rbind(x, (z[wa, ] + b2 * mz[wa] * mu) / com[wa] )
      n1 <- dim(x)[1]
    }
  }
  x[1:n, ]
}





