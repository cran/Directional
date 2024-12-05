################################
#### Hyper spherical-spherical regression
#### Tsagris Michail 11/2013
#### mtsagris@yahoo.gr
#### References: Chang Ted (1986)
#### Spherical egession. Annals of statistics, 14(3): 907-924
################################
hspher.reg <- function(y, x, xnew = NULL) {
  ## x is the independent variable
  ## y is the dependent variable
  d <- dim(y)[2]  ## dimensionality of the hyper-sphere
  xy <- crossprod(x, y)   ## crossprod(X, Y) / n, division by n not necessary though
  b <- svd(xy)  ## SVD of the XY matrix
  A <- tcrossprod(b$v, b$u )
  if ( det(A) < 0 ) {
    b$u[, d] <-  - b$u[, d]
    A <- tcrossprod(b$v, b$u )
  }

  est <- NULL
  if ( !is.null(xnew) )  est <- tcrossprod(xnew, A)

  list(A = A, est = est)
}
