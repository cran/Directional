dsespc <- function(y, mu, theta, logden = FALSE) {
  y <- as.matrix(y)
  if   ( dim(y)[2] == 1 )   y <- t(y)
  the1 <- theta[1]
  the2 <- theta[2]
  
  if ( the1 == 0  &  the2 == 0 ) {
    a <- as.vector( y %*% mu )
    rl <- sum(mu^2)
    d <- rl + 1 - a^2
    sqd <- sqrt(d)
    up <- log( (rl + 1) * sqd * ( atan2(sqd, -a) - atan2(sqd, a) + pi ) + 2 * a * d )
    down <- log(d^2)
    den <- up - down - log(4 * pi^2)
    if ( !logden )  den <- exp(den)

  } else { 
    heta <- sqrt(the1^2 + the2^2 + 1) - 1
    m0 <- sqrt(mu[2]^2 + mu[3]^2)
    rl <- sum(mu^2)  ## gamma^2
    x1b <- c( -m0^2, mu[1] * mu[2], mu[1] * mu[3] ) / ( m0 * sqrt(rl) )
    x2b <- c(0, -mu[3], mu[2])/m0
    T1 <- tcrossprod(x1b)
    T2 <- tcrossprod(x2b)
    T12 <- tcrossprod(x1b, x2b)
    vinv <- diag(3) + the1 * (T1 - T2) + the2 * ( T12 + t(T12) ) + 
            heta * (T1 + T2)
    a <- as.vector( y %*% mu )
    b <- Rfast::rowsums( y %*% vinv * y )
    E <- b * rl + b - a^2
    sqe <- sqrt(E)
    up <- log( b * (rl + 1) * sqe * ( atan2(sqe, -a) - atan2(sqe, a) + pi ) + 2 * a * E )
    down <- log( b * E^2 )
    den <- up - down -  log(4 * pi^2)
    if ( !logden )  den <- exp(den)
  }

  den
}
