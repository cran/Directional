iagd <- function(y, mu, logden = FALSE) {
  y <- as.matrix(y)
  p <- dim(y)[2]
  if ( p == 1 )   y <- t(y)
  a <- tcrossprod(mu, y)
  a2 <- a^2
  pa <- pnorm(a)
  da <- dnorm(a)
  gm <- pa + a2 * pa + a * da
  den <-  -log(2 * pi) + 0.5 * a2 - 0.5 * sum(mu^2) + log(gm)
  if ( !logden )  den <- exp(den)
  den
}
