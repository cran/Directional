dcircpurka <- function(x, m, a, rads = FALSE, logden = FALSE) {
  if ( !rads )  x <- x * pi/180
  x <- cbind( cos(x), sin(x) )
  m <- c( cos(m), sin(m) )
  den <- log(a) - log(2) - log(1 - exp(-a * pi)) - a * acos( x %*% m)
  if ( !logden )  den <- exp(den)
  den
}
