################################
#### Discrminant analysis for directional data
#### assuming a von Mises-Fisher distribution
#### Tsagris Michail 03/2014
#### mtsagris@yahoo.gr
#### References: J. E. Morris and P. J. Laycock (1974)
#### Discriminant Analysis of Directional Data (Biometrika)
################################
vmfda.pred <- function(xnew, x, ina) {
  ## xnew is the new observation(s)
  ## x is the data set
  ## ina is the group indicator variable
  xnew <- as.matrix(xnew)
  if ( ncol(xnew) == 1 )   xnew <- t(xnew)
  mod <- Rfast::multivmf.mle(x, ina, ell = FALSE)
  ki <- mod$ki
  p <- dim(x)[2]
  mat <- (p/2 - 1) * log(ki) + ki * tcrossprod(mod$mi, xnew) - log( besselI(ki, p/2 - 1, expon.scaled = TRUE) ) - ki
  Rfast::colMaxs(mat)
}
