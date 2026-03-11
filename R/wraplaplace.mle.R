wraplaplace.mle <- function(x, rads = FALSE) {
  if ( !rads )   x <- x * pi / 180
  n <- length(x)

  lik <- function(lam, x, n) {
    a <- ( exp( (2 * pi - x) * lam ) + exp(lam * x) ) / ( exp(2 * pi * lam) - 1 )
    n * log(lam) + sum( log(a) )
  }

  res <- optimize(lik, c(0, 150), x = x, n = n, maximum = TRUE)
  list( lambda = res$maximum, loglik = res$objective - n * log(2) )
}
