spcauchy.mle <- function(x, tol = 1e-06) {

  dm <- dim(x)
  n <- dm[1]  ;  d <- dm[2] - 1
  sp <- function(para, n, d) {
    rho <- 1 / ( 1 + exp( -para[1] ) )
    mu <- para[-1]
    mu <- mu / sqrt( sum(mu^2) )
    a <- as.vector(x %*% mu)
    - n * d * log(1 - rho^2) + d * sum( log1p( rho^2 - 2 * rho * a ) )
  }
  m1 <- optim( c( runif(1), rnorm(d + 1) ), sp, n = n, d = d, control = list(maxit = 5000) )
  m2 <- optim( c( runif(1), rnorm(d + 1) ), sp, n = n, d = d, control = list(maxit = 5000) )
  while (m1$value - m2$value > tol) {
    m1 <- m2
    m2 <- optim( m1$par, sp, n = n, d = d, control = list(maxit = 5000) )
  }
  rho <- 1 / ( 1 + exp( -m2$par[1] ) )
  mu <- m2$par[-1]
  mu <- mu / sqrt( sum(mu^2) )
  loglik <- n * lgamma( 0.5 * (d + 1) ) - 0.5 * n * (d + 1) * log(2 * pi)  - m2$value
  list(mu = mu, rho = rho, loglik = loglik)
}
