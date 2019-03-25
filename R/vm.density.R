vm.density <- function(x, k, m, rads = FALSE, logden = FALSE) {
  if ( !rads )  x <- x * pi/180
  den <- k * cos(x - m) - log(2 * pi) - log( besselI(k, 0, expon.scaled = TRUE) ) - k
  if ( !logden )  den <- exp(den)
  den
}
