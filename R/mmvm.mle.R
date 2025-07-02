mmvm.mle <- function(x, N, rads = FALSE) {
  if ( !rads )  x <- x * pi/180
  n <- length(x)

  likel <- function(par, x, N, n) {
    m <- par[1]  ;  k <- exp(par[2])
    - k * sum( cos(N * (x - m) ) ) + n * log(2 * pi) + n * ( log( besselI(k, 0, expon.scaled = TRUE) ) - k )
  }

  qa <- optim( c(mean(x), 1), likel, x = x, N = N, n = n, control = list(maxit = 5000) )
  list(mu = qa$par[1], kappa = exp(qa$par[2]), loglik = -qa$value )
}
