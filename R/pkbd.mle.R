pkbd.mle <- function(x, tol = 1e-6) {

  dm <- dim(x)
  n <- dm[1]  ;  d <- dm[2]

  pkbd <- function(rho, mu, x, n, d) {
    a <- as.vector(x %*% mu)
    n * log( 1 - rho^2 ) - 0.5 * d * sum( log1p( rho^2 - 2 * rho * a ) )
  }

  mu <- Rfast::colmeans(x)
  mu <- mu / sqrt(sum(mu^2) )
  mod <- optimize(pkbd, c(0, 1), mu = mu, x = x, n = n, d = d, maximum = TRUE, tol = 1e-6 )
  rho <- mod$maximum
  lik1 <- mod$objective

  down <- 1 + rho^2 - 2 * rho * as.vector( x %*% mu)
  mu <- Rfast::eachcol.apply(rho * x, down, oper = "/")
  mu <- mu / sqrt( sum(mu^2) )
  mod <- optimize(pkbd, c(0, 1), mu = mu, x = x, n = n, d = d, maximum = TRUE, tol = 1e-6 )
  rho <- mod$maximum
  lik2 <- mod$objective

  while ( abs( lik2 - lik1 ) > tol ) {
    lik1 <- lik2
    down <- 1 + rho^2 - 2 * rho * as.vector( x %*% mu)
    mu <- Rfast::eachcol.apply(rho * x, down, oper = "/")
    mu <- mu / sqrt( sum(mu^2) )
    mod <- optimize(pkbd, c(0, 1), mu = mu, x = x, n = n, d = d, maximum = TRUE, tol = 1e-6 )
    rho <- mod$maximum
    lik2 <- mod$objective
  }

  list(mu = mu, rho = rho, loglik = lik2 + n * lgamma(0.5 * d) - n * 0.5 * d * log(2 * pi) )
}
