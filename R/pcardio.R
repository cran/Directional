pcardio <- function(u, m, rho, rads = FALSE) {
  if ( !rads )  {
    u <- u * pi / 180
    m <- m * pi / 180
  }
  0.5 * u / pi + rho / pi * ( sin(u - m) + sin(m) )
}
