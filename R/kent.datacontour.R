################################
#### Contour plots of the Kent distribution on the sphere
#### with data appearing
#### Tsagris Michail 06/2014
#### mtsagris@yahoo.gr
################################
kent.datacontour <- function(x) {
  ## u contains the data in latitude and longitude
  ## the first column is the latitude and the
  ## second column is the longitude

  ## if u are eucliean coordinates turn them into
  ## latitude and longitude
  dm <- dim(x)
  if ( dm[2] == 3 ) {
    u <- euclid.inv(x)
  } else if ( dm[2] == 2 ) {
    u <- x
    x <- euclid(x) ## Euclidean coordinates used by Kent (1982)
  }

  n <- dm[1]  ## sample size
  a <- kent.mle(x) ## MLE estimation of the Kent distribution
  G <- a$G ## G matrix, the mean direction and the major-minor axes
  k <- a$para[1] ## kappa, concentration parameter
  b <- a$para[2] ## beta, ovalness
  gam <- c(0, k, 0)
  lam <- c(0, -b, b)
  ckb <- fb.saddle(gam, lam)[3] ## logarithm of the normalising constant
  n <- 100
  x1 <- seq(min(u[, 1]) - 5, max(u[, 1]) + 5, length = n)  ## latitude
  x2 <- seq(min(u[, 2]) - 5, max(u[, 2]) + 5, length = n)  ## longitude
  mat <- matrix(nrow = n, ncol = n)

  for (i in 1:n) {
    for (j in 1:n) {
      y <- euclid( c(x1[i], x2[j]) )
      can <-  -ckb + k * y %*% G[, 1] + b * (y %*% G[, 2])^2 -
        b * (y %*% G[, 3])^2
      if ( abs(exp( can) ) < Inf )  mat[i, j] <- exp(can)
    }
  }

  contour(x1, x2, mat, nlevels = 10, col = 2, xlab = "Latitude", ylab = "Longitude")
  points(u[, 1], u[, 2])
}
