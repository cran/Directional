pgcpc <- function(u, omega, g, rho, rads = FALSE) {
  if ( !rads )  {
    u <- u * pi / 180
    m <- m * pi / 180
  }

  funa <- function(u, omega, g, rho)  {
    uom <- u - omega
    a <- g * cos(uom)
    b <- cos(uom)^2 + sin(uom)^2 / rho
    0.5 / (pi * sqrt(rho) * ( b * sqrt(g^2 + 1) - a * sqrt(b) ) )
  }

  n <- length(u)
  lo <- rep(0, n)
  prob <- mapply( function(lo, u, omega, g, rho)  integrate(funa, lo, u, omega, g, rho)$value,
                                                  lo, u, omega = omega, g = g, rho = rho )
  prob
}
