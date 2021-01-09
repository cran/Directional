dpurka <- function(y, a, theta, logden = FALSE) {
  y <- as.matrix(y)
  p <- dim(y)[2]
  if ( p == 1 )   y <- t(y)
  p <- dim(y)[2]

  if ( p == 3 ) {
    A <- y %*% theta
    A[ abs(A) > 1 ] <- 1
    A <- acos(A)
    den <- log(a^2 + 1) - log(2 * pi) - log( 1 + exp( - a * pi ) ) - a * A
  } else if ( p > 3 ) {
    A <- y %*% theta
    A[ abs(A) > 1 ] <- 1
    A <- acos(A)
    den <- lgamma(p/2) - 0.5 * p * log(pi) + log(besselI(a, p - 1, expon.scaled = TRUE)) + a - a * A
  }

  if ( !logden )  den <- exp(den)
  den
}
