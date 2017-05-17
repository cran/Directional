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
  p <- dim(x)[2]  ## dimensonality of the data
  ina <- as.numeric(ina)
  g <- max(ina)
  mesi <- matrix(nrow = g, ncol = p)
  k <- numeric(g)
  nu <- nrow(xnew)
  mat <- matrix(nrow = nu, ncol = g)

  for (j in 1:g) {
    da <- Rfast::vmf.mle( x[ina == j, ] )  ## estimates the parameters of the vMF
    mesi[j, ] <- da$mu  ## mean direction
    k[j] <- da$kappa  ## concentration
    mat[, j] <- (p/2 - 1) * log(k[j]) + k[j] * xnew %*% mesi[j, ] - log( besselI(k[j], p/2 - 1, expon.scaled = TRUE) ) - k[j]   ##- 0.5 * p * log(2 * pi)
  }

  Rfast::rowMaxs(mat)
}
