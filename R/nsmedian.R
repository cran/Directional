nsmedian <- function(x, tol = 1e-07) {
  m <- Rfast::spat.med(x, tol = tol)
  m / sqrt( sum(m^2) )
}
