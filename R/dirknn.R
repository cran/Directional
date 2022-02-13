################################
#### k-NN algorithm for directional data
#### using the k-NN alorithm,
#### Tsagris Michail 01/2016
#### mtsagris@yahoo.gr
################################
dirknn <- function(xnew, ina, x, k = 5, type = "S", mesos = TRUE, parallel = FALSE, rann = FALSE) {
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
    if (rann) {
      for (m in 1:nc) {
        apo[[ m ]] <- Rfast::colSort( t( RANN::nn2(data = x[ina == m, ], query = xnew, k = max(k) )$nn.dists ) )
        apo[[ m ]] <-  0.5 * apo [[ m ]]^2 - 1
      }
    } else {
      for (m in 1:nc) {
        disa <- tcrossprod(x[ina == m, ], xnew)
        ## no need to compute the acos
        #disa[ disa >= 1 ] <- 1
        #disa[ disa <=  -1 ] <-  -1
        #adisa <- acos(disa)
        apo[[ m ]] <- Rfast::colSort(disa, descending = TRUE)[1:max(k), , drop = FALSE]
      }
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
    if ( rann ) {
      klen <- length(k)
      di <- RANN::nn2( data = x, query = xnew, k = max(k) )$nn.idx
      g <- matrix(nrow = nu, ncol = klen)
      m1 <- matrix(nrow = max(k), ncol = nu)
      for ( i in 1:nu )  m1[, i] <- ina[ di[i, ] ]
      for ( j in 1:klen ) g[, j] <- Rfast::colMaxs( Rfast::colTabulate( m1[1:k[j], ] ) )
    } else  g <- Rfast::dirknn(xnew, x, ina, k = k, type = "C", parallel = parallel)
  }  ## end if (type == "NS")
  colnames(g) <- paste("k=", k, sep = "")
  g
}
