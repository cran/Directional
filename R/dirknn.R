################################
#### Discriminant analysis for directional data
#### using the k-NN alorithm,
#### Tsagris Michail 01/2016
#### mtsagris@yahoo.gr
################################
dirknn <- function(x, xnew, k = 5, ina, type = "S", mesos = TRUE) {
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


  if (type == "NS") {
    ## Non Standard algorithm
    klen <- length(k)
    g <- matrix(0, nu, klen)
    ta <- matrix(nrow = nu, ncol = nc)
    apo <- list()
    for (m in 1:nc) {
	    disa <- tcrossprod(x[ina == m, ], xnew)
      disa[ disa >= 1 ] <- 1
      disa[ disa <=  -1 ] <-  -1
      disa <- acos(disa)
      apo[[ m ]] <- Rfast::sort_mat(disa)[1:max(k), , drop = FALSE]
    }
    for (j in 1:klen) {
      for (m in 1:nc) {
        if ( mesos ) {
          ta[, m] <- Rfast::colmeans( apo[[ m ]][1:k[j], , drop = FALSE] )
        } else  ta[, m] <- Rfast::colhameans( apo[[ m ]][1:k[j], , drop = FALSE] )
      }
      g[, j] <- Rfast::rowMins(ta)
    }

  } else {   ## if type is "S"   ## Standard algorithm
    g <- Rfast::dirknn(xnew, x, ina, k = k, type = "C", parallel = FALSE)
  }  ## end if (type == "NS")
  colnames(g) <- paste("k=", k, sep = "")
  g
}
