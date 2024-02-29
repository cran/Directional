pwrapcauchy <- function(u, m, rho, rads = FALSE) {
  if ( !rads )  {
    u <- u * pi / 180
    m <- m * pi / 180
  }

  f <- 0.5 * (1 - rho^2) / pi
  funa <- function(u, m, rho, f)    f / (1 + rho^2 - 2 * rho * cos(u - m) )
  n <- length(u)
  lo <- rep(0, n)
  prob <- mapply( function(lo, u, m, rho, f)  integrate(funa, lo, u, m, rho, f)$value, lo, u, m = m, rho = rho, f = f )
  prob
}
