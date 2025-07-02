circbeta.mle <- function(x, rads = FALSE) {
  if ( !rads )   x <- x * pi/180
  n <- length(x)

  fun <- function(par, x) {
    m <- cos(par[1])
    a <- exp(par[2])
    b <- exp(par[3])
    con <- cos(x - m)
    den <-  -(a + b) * log(2) - lbeta(a, b) + (a - 0.5) * log1p(con) + (b - 0.5) * log(1 - con)
    -sum(den)
  }

  qa <- optim( c( mean(x), rnorm(2) ), fun, x = x, control = list(maxit = 5000) )
  qa <- optim( qa$par, rnorm(2), fun, x = x, control = list(maxit = 5000) )
  param <- c( cos(qa$par[1]), exp(qa$par[2]), exp(qa$par[3]) )
  names(param) <- c("mean", "alpha", "beta")
  list( loglik = -qa$value, param = param)
}
