pspml <- function(u, mu, rads = FALSE) {
  if ( !rads )  u <- u * pi / 180
  gam <- sqrt( sum(mu^2) )
  m <- ( atan(mu[2] / mu[1]) + pi * I(mu[1] < 0) ) %% (2 * pi)
  f <- 0.5 * exp(-0.5 * gam^2) / pi

  funa <- function(u, f, gam) {
    com <- gam * cos(u - m)
    f * (1 + com * pnorm(com) / dnorm(com) )
  }
  n <- length(u)
  lo <- rep(0, n)
  prob <- mapply( function(lo, u, f, gam)  integrate(funa, lo, u, f, gam)$value, lo, u, f = f, gam = gam )
  prob
}
