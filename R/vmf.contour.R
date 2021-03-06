################################
#### Contour plots of the von Mises-Fisher distribution on the sphere
#### Tsagris Michail 06/2014
#### mtsagris@yahoo.gr
################################

vmf.contour <- function(k) {
  ## k is the concentration parameter
  rho <- pi/2  ## radius of the circular disc
  x <- seq(-rho, rho, by = 0.01)
  n <- length(x)
  mat <- matrix(rep(x^2, n), ncol = n)
  z <- mat + t(mat)
  theta <- sqrt(z)
  ind <- ( theta < rho )  ## checks if x^2+y^2 < rho^2
  xa <- 0.5 * log(k) + k * cos(theta) - 1.5 * log(2 * pi) - log( besselI(k, 0.5, expon.scaled = TRUE) ) - k
  mat <- exp(xa) * ind
  contour(x, x, mat)
}
