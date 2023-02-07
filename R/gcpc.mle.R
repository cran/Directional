gcpc.mle <- function(x, rads = FALSE) {

  lik <- function(param, x, n) {
    mu <- param[1:2]
    rho <- 1 / ( 1 + exp(-param[3]) )
    g2 <- sum(mu^2)
    ksi <- mu / sqrt(g2)
    sinv <- matrix( c( ksi[1]^2 + ksi[2]^2/rho, ksi[1] * ksi[2] * (1 - 1/rho),
                       ksi[1] * ksi[2] * (1 - 1/rho), ksi[2]^2 + ksi[1]^2/rho ), ncol = 2 )
    a <- as.vector(x %*% mu)
    b <- rowsums( x %*% sinv * x )
    n * 0.5 * log(rho) + sum( log( b * sqrt(g2 + 1) - a * sqrt(b) ) )
  }

  if ( !rads )  x <- x * pi/180
  x <- cbind( cos(x), sin(x) )
  n <- dim(x)[1]

  mod <- optim( rnorm(3), lik, x = x, n = n, control = list(maxit = 5000) )
  lik1 <-  -mod$value
  mod <- optim( mod$par, lik, x = x, n = n, control = list(maxit = 5000) )
  lik2 <-  -mod$value
  while (lik2 - lik1 > 1e-06) {
    lik1 <- lik2
    mod <- optim( mod$par, lik, x = x, n = n, control = list(maxit = 5000) )
  }

  mu <- mod$par[1:2]
  rho <- 1 / ( 1 + exp(-mod$par[3]) )
  gama <- sqrt( sum(mu^2) )
  circmu <- ( atan(mu[2]/mu[1]) + pi * I(mu[1] < 0) ) %% (2 * pi)
  list(mu = mu, circmu = circmu, gamma = gama, rho = rho, loglik = -mod$value - n * log(2 * pi) )
}
