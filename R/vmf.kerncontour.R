################################
#### Contour plots of the von Mises-Fisher kernel density estimate on the sphere
#### Tsagris Michail 03/2015
#### mtsagris@yahoo.gr
#### Garcia-Portugues E. (2013)
#### Exact risk improvement of bandwidth selectors for kernel
#### density estimation with directional data
#### Electronic Journal of Statistics
################################

vmf.kerncontour <- function(u, thumb = "none", den.ret = FALSE, full = FALSE,
                            ngrid = 100) {
  ## u contains the data in latitude and longitude
  ## the first column is the latitude and the
  ## second column is the longitude
  ## thumb is either 'none' (defualt), or 'rot' (Garcia-Portugues, 2013)
  ## den.ret if set to TRUE returns a list with the following components:
  ##  * lat - latitudes of densities
  ##  * long - longitudes of densities
  ##  * h - bandwidth used in calculation
  ##  * den - matrix with densities at each latitude / longitude
  ## full if set to TRUE calculates densities for the full sphere, otherwise
  ##   using extents of the data
  ## ngrid specifies the number of points taken at each axis
  n <- dim(u)[1]  ## sample size
  x <- euclid(u)

  if (thumb == "none") {
    h <- as.numeric( vmfkde.tune(x, low = 0.1, up = 1)[1] )
  } else if (thumb == "rot") {
    k <- vmf.mle(x)$kappa
    h <- ( (8 * sinh(k)^2) / (k * n * ( (1 + 4 * k^2) * sinh(2 * k) -
    2 * k * cosh(2 * k)) ) ) ^ ( 1/6 )
  }

  if (full) {
    x1 <- seq( 0, 180, length = ngrid )  ## latitude
    x2 <- seq( 0, 360, length = ngrid )  ## longitude
  } else {
    x1 <- seq( min(u[, 1]) - 5, max(u[, 1]) + 5, length = ngrid )  ## latitude
    x2 <- seq( min(u[, 2]) - 5, max(u[, 2]) + 5, length = ngrid )  ## longitude
  }
  cpk <- 1 / ( ( h^2)^0.5 *(2 * pi)^1.5 * besselI(1/h^2, 0.5) )
  mat <- matrix(nrow = ngrid, ncol = ngrid)

  for (i in 1:ngrid) {
    for (j in 1:ngrid) {
      y <- Directional::euclid( c(x1[i], x2[j]) )
      a <- as.vector( tcrossprod(x, y / h^2) )
      can <- sum( exp(a + log(cpk)) ) / ngrid
      if (abs(can) < Inf)   mat[i, j] <- can
    }
  }

  if (den.ret) {
    return(list(lat = x1, long = x2, h = h, den = mat))
  } else {
    contour(x1, x2, mat, nlevels = 10, col = 2, xlab = "Latitude",
            ylab = "Longitude")
    points(u[, 1], u[, 2])
  }
}
