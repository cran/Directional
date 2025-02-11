wrapnormal.mle <- function(x, rads = FALSE) {
  if ( !rads )   x <- x * pi/180
  y <- cbind( cos(x), sin(x) )
  n <- dim(y)[1]
  m <- Rfast::colmeans(y)
  rho <- sqrt( sum(m^2) )
  mu <- ( atan(m[2]/m[1]) + pi * I(m[1] < 0) ) %% (2 * pi)
  y <- x - mu
  k <- 1:100
  y <- Rfast::rep_row(y, 100)
  y <- Rfast::rowsums( rho^(k^2) * cos(k * y) )
  loglik <-  - n * log(2 * pi) + sum( log1p(2 * y) )
  param <- c(mu, rho)
  names(param) <- c("mean direction", "rho")
  list(loglik = loglik, param = param)
}
