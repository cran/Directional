dwrapnormal <- function(x, m, rho, rads = FALSE, logden = FALSE) {
  if ( !rads )   x <- x * pi/180
  n <- length(x)
  y <- x - m
  k <- 1:100
  y <- Rfast::rep_row(y, 100)
  y <- Rfast::rowsums( rho^(k^2) * cos(k * y) )
  den <-  -log(2 * pi) + log1p(2 * y)
  if ( !logden )  den <- exp(den)
  den
}
