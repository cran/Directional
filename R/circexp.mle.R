circexp.mle <- function(x, rads = FALSE, tol = 1e-06) {
  if ( !rads )   x <- x * pi/180
  n <- length(x)
  sx <- sum(x)
  fun <- function(lam, n, sx)  n * log(lam) - lam * sx - sum( log( 1 - exp(-2 * pi * lam) ) )
  mod <- optimise(fun, c(0.00001, 1000), n = n, sx = sx, tol = tol, maximum = TRUE)
  list(loglik = mod$objective, lambda = mod$maximum)
}
