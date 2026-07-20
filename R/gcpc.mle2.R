gcpc.mle2 <- function(x, rads = FALSE) {
  lik <- function(para, x, n) {
    omega <- para[1]
    r <- exp(para[2])  ;  g <- exp(para[3])
    phi <- x - omega
    a <- g * cos(phi)
    b <- cos(phi)^2 + sin(phi)^2/r
    0.5 * n * log(r) + sum( log( b * sqrt(g^2 + 1) - a * sqrt(b) ) )
  }

  if ( !rads )  x <- x / 180 * pi
  n <- length(x)
  mod <- Directional::cipc.mle(x, TRUE)
  ini <- c(mod$circmu, log( mod$gamma), rnorm(1) )
  suppressWarnings({
    mod <- optim(ini, lik, x = x, n = n)
    mod <- optim(mod$par, lik, x = x, n = n)
    mod <- optim(mod$par, lik, x = x, n = n)
  })
  circmu <- mod$par[1] %% (2 *pi)
  rho <- exp(mod$par[2])  ;  gama <- exp(mod$par[3])
  mu <- c( cos(circmu), sin(circmu) ) * gama
  list(mu = mu, circmu = circmu, gamma = gama, rho = rho, loglik = -mod$value - n * log(2 * pi) )
}
