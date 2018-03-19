################################
#### Spherical-spherical regression
#### Tsagris Michail 11/2013
#### mtsagris@yahoo.gr
#### References: Chang Ted (1986)
#### Spherical egession. Annals of statistics, 14(3): 907-924
################################
spher.reg <- function(y, x, rads = FALSE) {
  ## x is the independent variable
  ## y is the dependent variable
  ## The first row of both matrices is the latitude
  ## and the second is the longitude
  n <- dim(x)[1]  ## sample size
  if ( dim(x)[2] == 2  &  dim(y)[2] == 2 ) {
    if ( !rads ) {
      x <- pi * x / 180  ## from degrees to rads
      y <- pi * y / 180
    }  ## from degrees to rads
    ## the first row of both matrices is the latitude and the second is the longitude
    ## the next two rows transform the data to Euclidean coordinates
    cosx1 <- cos(x[, 1])   ;  cosy1 <- cos(y[, 1])
    X <- cbind( cosx1 * cos(x[, 2]), cosx1 * sin(x[, 2]), sin(x[, 1]) )
    Y <- cbind( cosy1 * cos(y[, 2]), cosy1 * sin(y[, 2]), sin(y[, 1]) )
  } else if ( dim(x)[2] == 3  &  dim(y)[2] == 3 ) {
    X <- x
    Y <- y
  }

  XY <- crossprod(X, Y) / n
  b <- svd(XY)  ## SVD of the XY matrix
  A <- tcrossprod(b$v, b$u )
  if ( det(A) < 0 ) {
    b$u[, 3] <-  - b$u[, 3]
    A <- tcrossprod(b$v, b$u )
  }

  est <- tcrossprod(X, A)
  list(A = A, fitted = est)
}
