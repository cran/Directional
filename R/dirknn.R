################################
#### k-NN algorithm for directional data
#### using the k-NN alorithm,
#### Tsagris Michail 01/2016
#### mtsagris@yahoo.gr
################################
dirknn <- function(xnew, ina, x, k = 5, mesos = TRUE, parallel = FALSE, rann = FALSE) {
  ## x is the matrix containing the data
  ## xnew is the new data
  ## k is the number of neighbours to use
  ## ina indicates the groups, numerical variable
  ## type is either 'S' or 'NS'. Should the standard k-NN be use or not
  ## if mesos is TRUE, then the arithmetic mean distange of the k nearest
  ## points will be used.
  ## If not, then the harmonic mean will be used. Both of these apply for
  ## the non-standard algorithm, that is when type='NS'
  xnew <- matrix(xnew, ncol = dim(x)[2])  ## makes sure xnew is a matrix
  ina <- as.numeric(ina) ## makes sure ina is numeric
  nc <- max(ina)  ## The number of groups
  nu <- dim(xnew)[1]

  if ( rann ) {
    klen <- length(k)
    di <- Rnanoflann::nn( data = x, points = xnew, k = max(k) )$indices
    g <- matrix(nrow = nu, ncol = klen)
    m1 <- matrix(nrow = max(k), ncol = nu)
    for ( i in 1:nu )  m1[, i] <- ina[ di[i, ] ]
    for ( j in 1:klen ) g[, j] <- Rfast::colMaxs( Rfast::colTabulate( m1[1:k[j], ] ) )
  } else  g <- Rfast::dirknn(xnew, x, ina, k = k, type = "C", parallel = parallel)
  colnames(g) <- paste("k=", k, sep = "")
  g
}
