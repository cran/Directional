dcipc <- function(x, omega, g, rads = FALSE, logden = FALSE) {
  if (!rads)  x <- x * pi/180
  den <-  - log(2 * pi) - log( sqrt(g^2 + 1) - g * cos(x - omega) )
  if ( !logden )  den <- exp(den)
  den
}
