cipc.mle <- function(x, rads = FALSE) {

  lik <- function(mu, x) {
    g2 <- sum(mu^2)
    a <- as.vector(x %*% mu)
    sum( log( sqrt(g2 + 1) - a ) )
  }

  if ( !rads )  x <- x * pi/180
  x <- cbind( cos(x), sin(x) )
  n <- dim(x)[1]

  mod <- optim( rnorm(2), lik, x = x, control = list(maxit = 5000) )
  mod <- optim( mod$par, lik, x = x, control = list(maxit = 5000) )
  mu <- mod$par
  circmu <- ( atan(mu[2]/mu[1]) + pi * I(mu[1] < 0) ) %% (2 * pi)
  gama <- sqrt( sum(mu^2) )
  list(mu = mu, circmu = circmu, gamma = gama, loglik = -mod$value - n * log(2 * pi) )
}
