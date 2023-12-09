################################
#### Discriminant analysis for directional data
#### using the k-NN alorithm,
#### Tsagris Michail 01/2016
#### mtsagris@yahoo.gr
################################
cosnn <- function(xnew, x, k = 5, index = FALSE, rann = FALSE) {
  xnew <- matrix(xnew, ncol = dim(x)[2])  ## makes sure xnew is a matrix
  nu <- dim(xnew)[1]
  if ( index ) {
    if (rann) {
      disa <-  t( Rnanoflann::nn( data = x, points = xnew, k = k )$indices )
    } else disa <- Rfast::dista(xnew, x, k = k, square = TRUE, index = TRUE )
  } else {
    if (rann) {
      disa <-  t( Rnanoflann::nn(data = x, points = xnew, k = k, square = TRUE )$distances )
    } else disa <- Rfast::dista(xnew, x, k = k, square = TRUE)
    disa <- Rfast::colSort( 0.5 * disa - 1)
    disa[ disa >= 1 ] <- 1
    disa[ disa <=  -1 ] <-  -1
    disa <- acos(disa)
  }
  disa
}
