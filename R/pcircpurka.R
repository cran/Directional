pcircpurka <- function(u, m, a, rads = FALSE) {
  if ( !rads )  {
    u <- u * pi / 180
    m <- m * pi / 180
  }
  a1 <- (u - m) * ( exp(a * abs(u - m) ) - 1 ) * exp(a * pi - a * abs(u - m) )
  b1 <- 2 * ( exp(a * pi) - 1 ) * abs(u - m)
  a0 <- - m * ( exp(a * m) - 1 ) * exp(a * pi - a * m )
  b0 <- 2 * ( exp(a * pi) - 1 ) * m
  a1 / b1 - a0/b0
}

