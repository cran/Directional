pcircbeta <- function(u, m, a, b, rads = FALSE) {
  if ( !rads )  {
    u <- u * pi / 180
    m <- m * pi / 180
  }

  down <- 1 / ( 2^(a + b) * beta(a, b) )
  funa <- function(u, m, a, b, down) {
    com <- cos(u - m)
    down * ( 1 + com )^(a - 0.5) * ( 1 - com )^(b - 0.5)
  }
  n <- length(u)
  lo <- rep(0, n)
  prob <- mapply( function(lo, u, m, a, b, down)  integrate(funa, lo, u, m, a, b, down)$value,
                                                  lo, u, m = m, a = a, b = b, down = down )
  prob
}
