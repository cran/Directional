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
  x <- as.matrix(x)  ## makes sure the x is a matrix
  x <- x / sqrt( rowSums(x^2) )  ## makes sure x are unit vectors
  xnew <- matrix(xnew, ncol = ncol(x))  ## makes sure xnew is a matrix
  xnew <- xnew / sqrt( rowSums(xnew^2) ) ## makes sure xnew are unit vectors
  n <- nrow(x)  ## sample size
  ina <- as.numeric(ina) ## makes sure ina is numeric
  nc <- max(ina)  ## The number of groups
  nu <- nrow(xnew)
  z <- rbind(x, xnew)
  nz <- nrow(z)
  dis <- crossprod( t(z) )
  dis <- as.matrix(dis)
  diag(dis) <- 1
  dis <- acos(dis)
  diag(dis) <- 0 
  g <- rep(0, nu)
  ta <- matrix(nrow = nu, ncol = nc)
  apo <- matrix( dis[-c(1:n), 1:n], nrow = nu )
  if (type == "NS") {
    ## Non Standard algorithm
    for (l in 1:nu) {
      for (m in 1:nc) {
        dista <- apo[l, ina == m]
        if (mesos == TRUE) {
          ta[l, m] <- mean( sort(dista)[1:k] )
        } else {
          ta[l, m] <- k / sum( 1 / sort(dista)[1:k] )
        }
      }
    }
    g <- apply(ta, 1, which.min)
  }
  if (type == "S") {
    ## Standard algorithm
    for (l in 1:nu) {
      xa <- cbind(ina, apo[l, ])
      qan <- xa[order(xa[, 2]), ]
      sa <- qan[1:k, 1]
      tab <- table(sa)
      g[l] <- as.integer( names(tab)[ which.max(tab) ] )
    }
  }
  return(g)
}