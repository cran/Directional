dcircbeta <- function(x, m, a, b, rads = FALSE, logden = FALSE) {
  if ( !rads )  x <- x * pi/180
  con <- cos(x - m)
  den <-  -(a + b) * log(2) - lbeta(a, b) + (a - 0.5) * log(1 + con) + (b - 0.5) * log(1 - con)
  if ( !logden )  den <- exp(den)
  den
}
