gcpc.mle <- function(x, rads = FALSE) {

  likint <- function(mu, rho, x, n) {
    g2 <- sum(mu^2)
    ksi <- mu / sqrt(g2)
    sinv <- matrix( c( ksi[1]^2 + ksi[2]^2/rho, ksi[1] * ksi[2] * (1 - 1/rho),
                       ksi[1] * ksi[2] * (1 - 1/rho), ksi[2]^2 + ksi[1]^2/rho ), ncol = 2 )
    a <- as.vector(x %*% mu)
    b <- Rfast::rowsums( x %*% sinv * x )
    n * 0.5 * log(rho) + sum( log( b * sqrt(g2 + 1) - a * sqrt(b) ) )
  }

  lik0 <- function(rho, x, n, ma) {
    m1 <- optim(ma, likint, rho = rho, x = x, n = n, control = list(maxit = 1000) )
    - optim(m1$par, likint, rho = rho, x = x, n = n, control = list(maxit = 1000) )$value
  }

  lik <- function(param, x, n) {
    mu <- param[1:2]
    rho <- param[3]
    g2 <- sum(mu^2)
    ksi <- mu / sqrt(g2)
    sinv <- matrix( c( ksi[1]^2 + ksi[2]^2/rho, ksi[1] * ksi[2] * (1 - 1/rho),
                       ksi[1] * ksi[2] * (1 - 1/rho), ksi[2]^2 + ksi[1]^2/rho ), ncol = 2 )
    a <- as.vector(x %*% mu)
    b <- Rfast::rowsums( x %*% sinv * x )
    n * 0.5 * log(rho) + sum( log( b * sqrt(g2 + 1) - a * sqrt(b) ) )
  }

  ma <- Directional::cipc.mle(x, rads = rads)$mu
  if ( !rads )  x <- x * pi/180
  x <- cbind( cos(x), sin(x) )
  n <- dim(x)[1]

  rho <- optimise(lik0, c(0.001, 1000), x = x, n = n, ma = ma, maximum = TRUE)$maximum
  mod <- optim(ma, likint, rho = rho, x = x, n = n, control = list(maxit = 5000) )
  suppressWarnings({
    mod <- optim(c(mod$par, rho), lik, x = x, n = n, control = list(maxit = 5000) )
  })
  mu <- mod$par[1:2]  ;  rho = mod$par[3]
  gama <- sqrt( sum(mu^2) )
  circmu <- ( atan(mu[2]/mu[1]) + pi * I(mu[1] < 0) ) %% (2 * pi)
  list(mu = mu, circmu = circmu, gamma = gama, rho = rho, loglik = -mod$value - n * log(2 * pi) )
}
