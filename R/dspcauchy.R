dspcauchy <- function(y, mu, rho, logden = FALSE) {

  y <- as.matrix(y)
  if ( dim(y)[2] == 1 )  y <- t(y)
  d <- dim(y)[2] - 1

  a <- as.vector(y %*% mu)
  den <-  d * log(1 - rho^2) - d * log1p( rho^2 - 2 * rho * a ) +
          lgamma(0.5 * (d + 1)) - 0.5 * (d + 1) * log(pi) - log(2)
  if ( !logden )  den <- exp(den)
  den
}
