rcircbeta <- function(n, m, a, b, rads = TRUE) {
  x <- rangen::Rbeta(n, a, b)
  y <- acos( 2 * x - 1 )
  u <- rangen::Runif(n)
  x <- y * (u < 0.5) + (2 * pi - y) * (u > 0.5) + m
  if (!rads)  x <- x/pi * 180
  x
}
