sespc.mle <- function(y, full = FALSE, tol = 1e-6) {

  n <- dim(y)[1]
  I3 <- diag(3)
   mag <- function(param, y) {
     m <- param[1:3]
     the1 <- param[4]
     the2 <- param[5]
     heta <- sqrt(the1^2 + the2^2 + 1) - 1
     m0 <- sqrt(m[2]^2 + m[3]^2)
     rl <- sum(m^2)  ## gamma^2
     x1b <- c( -m0^2, m[1] * m[2], m[1] * m[3] ) / ( m0 * sqrt(rl) )
     x2b <- c(0, -m[3], m[2])/m0
     T1 <- tcrossprod(x1b)
     T2 <- tcrossprod(x2b)
     T12 <- tcrossprod(x1b, x2b)
     vinv <- I3 + the1 * (T1 - T2) + the2 * ( T12 + t(T12) ) +
             heta * (T1 + T2)
     a <- as.vector( y %*% m )
     b <- rowSums( y %*% vinv * y )
     E <- b * rl + b - a^2
     sqe <- sqrt(E)
     up <- log( b * (rl + 1) * sqe * ( atan2(sqe, -a) - atan2(sqe, a) + pi ) + 2 * a * E )
     down <- log( b * E^2 )
     - sum(up) + sum(down)
   }
  mod <- Directional::sipc.mle(y)
  da <- nlm( mag, c( mod$mu, rnorm(2) ), y = y, iterlim = 10000 )
  lik1 <-  -da$minimum
  da <- optim( da$estimate, mag, y = y, control = list(maxit = 10000) )
  lik2 <-  -da$value
  while (lik2 - lik1 > tol) {
    lik1 <- lik2
    da <- optim( da$par, mag, y = y, control = list(maxit = 10000) )
    lik2 <-  -da$value
  }

  if (full) {
    m <- da$par[1:3]
    the1 <- da$par[4]
    the2 <- da$par[5]
    theta <- sqrt(the1^2 + the2^2 + 1) - 1
    m0 <- sqrt(m[2]^2 + m[3]^2)
    rl <- sum(m^2)  ## gamma^2
    x1b <- c( -m0^2, m[1] * m[2], m[1] * m[3] ) / ( m0 * sqrt(rl) )
    x2b <- c(0, -m[3], m[2])/m0
    T1 <- tcrossprod(x1b)
    T2 <- tcrossprod(x2b)
    T12 <- tcrossprod(x1b, x2b)
    vinv <- I3 + the1 * (T1 - T2) + the2 * ( T12 + t(T12) ) + theta * (T1 + T2)
    lam2 <- theta + 1 - 0.5 * sqrt( (2 * theta + 2 )^ 2 - 4 )
    psi <- 0.5 * acos( 2 * the1 / (1/lam2 - lam2 ) )
    res <- list(mu = m, theta = c(the1, the2), loglik = lik2 - n * log(4 * pi^2),
                vinv = vinv, lambda = lam2, psi = psi, sipc.loglik = mod$loglik)
  } else  res <- list(mu = da$par[1:3], theta = da$par[4:5], loglik = lik2 - n * log(4 * pi^2), sipc.loglik = mod$loglik )

  res
}
