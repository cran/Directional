pcipc <- function(u, omega, g, rads = FALSE) {
  if ( !rads )  {
    u <- u * pi / 180
    m <- m * pi / 180
  }

  funa <- function(u, omega, g)  0.5 / ( pi * ( sqrt(g^2 + 1) - g * cos(u - omega) ) )

  n <- length(u)
  lo <- rep(0, n)
  prob <- mapply( function(lo, u, omega, g)  integrate(funa, lo, u, omega, g)$value, lo, u, omega = omega, g = g )
  prob
}
