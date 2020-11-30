rspml <- function(n, mu, rads = TRUE) {
  y <- Rfast::matrnorm(n, 3) + rep(mu, rep(n, 2))
  y / sqrt( Rfast::rowsums(y^2) )  ## ESAG simulated values
  y <- ( atan(y[, 2]/y[, 1]) + pi * I(y[, 1] < 0) ) %% (2 * pi)
  if ( !rads )  y <- y * pi/180
  y
}
