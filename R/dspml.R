dspml <- function(x, mu, rads = FALSE, logden = FALSE) {
  if ( !rads )  x <- x * pi/180
  x <- cbind( cos(x), sin(x) )
  ta <- tcrossprod(mu, x)
  den <-  -0.5 * sum(mu^2) + log1p( ta * pnorm(ta)/dnorm(ta) ) - log(2 * pi)
  if ( !logden )  den <- exp(den)
  den
}
