pmmvm <- function(u, m, k, N, rads = FALSE) {
  if ( !rads )  {
    u <- u * pi / 180
    m <- m * pi / 180
  }
  if ( k > 0 ) {
    f <- 2 * pi * besselI(k, 0)
    funa <- function(u, m, k, N)  exp(k * cos(N * (u - m) ) )
    n <- length(u)
    lo <- rep(0, n)
    prob <- mapply( function(lo, u, m, k, N)  integrate(funa, lo, u, m, k, N)$value, lo, u, m = m, k = k, N = N) / f
  } else  prob <- u / ( 2 * pi )
  prob
}


