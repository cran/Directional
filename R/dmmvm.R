dmmvm <- function(x, m, k, N, rads = FALSE, logden = FALSE) {
  if ( !rads )  x <- x * pi/180
  den <- k * cos(N * (x - m) ) - log(2 * pi) - log( besselI(k, 0, expon.scaled = TRUE) ) - k
  if ( !logden )  den <- exp(den)
  den
}
