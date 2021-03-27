dspml <- function(x, mu, rads = FALSE, logden = FALSE) {
  if ( !rads )  x <- x * pi/180
  x <- cbind( cos(x), sin(x) )
  tau <- tcrossprod(mu, x)
  ptau <- pnorm(tau) 
  gam <- sum(mu^2)
  f <-  - 0.5   ;   con <- sqrt(2 * pi) 
  den <-  - 0.5 * gam + log1p( tau * ptau * con / exp(f * tau^2) ) - log(2 * pi)
  if ( !logden )  den <- exp(den)
  den
}
