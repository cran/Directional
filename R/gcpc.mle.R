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

  lik0 <- function(rho, x, n) {
    m1 <- optim(rnorm(2), likint, rho = rho, x = x, n = n, control = list(maxit = 1000) )
    - optim(m1$par, likint, rho = rho, x = x, n = n, control = list(maxit = 1000) )$value
  }

 # lik <- function(param, x, n) {
 #   mu <- param[1:2]
 #   rho <- 1 / ( 1 + exp(-param[3]) )
 #   g2 <- sum(mu^2)
 #   ksi <- mu / sqrt(g2)
 #   sinv <- matrix( c( ksi[1]^2 + ksi[2]^2/rho, ksi[1] * ksi[2] * (1 - 1/rho),
 #                      ksi[1] * ksi[2] * (1 - 1/rho), ksi[2]^2 + ksi[1]^2/rho ), ncol = 2 )
 #   a <- as.vector(x %*% mu)
 #   b <- Rfast::rowsums( x %*% sinv * x )
 #   n * 0.5 * log(rho) + sum( log( b * sqrt(g2 + 1) - a * sqrt(b) ) )
 # }

  if ( !rads )  x <- x * pi/180
  x <- cbind( cos(x), sin(x) )
  n <- dim(x)[1]

  rho <- optimise(lik0, c(0.001, 1), x = x, n = n, maximum = TRUE)$maximum
  mod <- optim(rnorm(2), likint, rho = rho, x = x, n = n, control = list(maxit = 5000) )
  mod <- optim(mod$par, likint, rho = rho, x = x, n = n, control = list(maxit = 5000) )
  mu <- mod$par
  gama <- sqrt( sum(mu^2) )
  circmu <- ( atan(mu[2]/mu[1]) + pi * I(mu[1] < 0) ) %% (2 * pi)
  list(mu = mu, circmu = circmu, gamma = gama, rho = rho, loglik = -mod$value - n * log(2 * pi) )
}
