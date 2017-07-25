wrapcauchy <- function(x, rads = FALSE, tol = 1e-07) {
  if ( !rads )   x <- x * pi/180
  Rfast::wrapcauchy.mle(x, tol)
}








