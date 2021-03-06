################################
#### Contour plots of the Kent distribution on the sphere
#### Tsagris Michail 06/2014
#### mtsagris@yahoo.gr
################################

kent.contour <- function(k, b) {
  ## k is the concentration parameter
  ## b is the ovalness parameter
  ## b must be less than k/2
  gam <- c(0, k, 0)
  lam <- c(0, -b, b)
  con <- Directional::fb.saddle(gam, lam)[3]
  rho <- sqrt(2)
  x <- seq(-rho, rho, by = 0.01)
  n <- length(x)
  mat1 <- matrix(rep(x^2, n), ncol = n)
  mat2 <- t(mat1)
  z <- sqrt( mat1 + mat2 )
  ind <- ( z^2 < rho^2 )  ## checks if x^2+y^2 < rho^2
  ind[ !ind ] <- NA
  theta <- 2 * asin(0.5 * z)
  xa <- k * cos(theta) + b * (mat1 - mat2) - con
  mat <- exp(xa) * ind
  contour(x, x, mat)
}
