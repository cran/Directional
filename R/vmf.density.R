vmf.density <- function(y, k, mu, logden = FALSE) {
  y <- as.matrix(y)
  p <- dim(y)[2]
  if ( p == 1 )   y <- t(y)
  p <- dim(y)[2]
  den <- (p/2 - 1) * log(k) - 0.5 * p * log(2 * pi) + k * tcrossprod(mu, y) - log( besselI(k, p/2 - 1, expon.scaled = TRUE) ) - k
  if ( !logden )  den <- exp(den)
  den
}
