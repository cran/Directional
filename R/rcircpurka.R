rcircpurka <- function(n, m, a, rads = TRUE) {

  u <- runif(n)
  con <- ( 2 - exp(-a * m) - exp(a * m - 2 * pi * a) ) / a
  u1 <- u[u < 0.5]
  x1 <- log(1 + a / con * u1 * exp(a * m) ) / a
  u2 <- u[u > 0.5]
  x2 <- m - log(1 - a / con * (u2 - 0.5) ) / a
  x <- c(x1, x2)
  if ( !rads )  x <- x/pi * 180
  x

}
