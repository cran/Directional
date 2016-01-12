################################
#### Discriminant analysis for directional data
#### using the k-NN alorithm, tuning the k neighbours
#### Tsagris Michail 01/2016 
#### mtsagris@yahoo.gr
################################

dirknn.tune <- function(x, M = 10, A = 5, ina, type = "S", 
  mesos = TRUE, seed = TRUE) {
  ## x is the matrix containing the data
  ## M is the number of folds, set to 10 by default
  ## A is the maximum number of neighbours to use
  ## actually it goes until A + 1 
  ## ina indicates the groups, numerical variable
  ## type is either 'S' or 'NS'. Should the standard k-NN be use or not
  ## if mesos is TRUE, then the arithmetic mean distange of the k nearest 
  ## points will be used.
  ## If not, then the harmonic mean will be used. Both of these apply for 
  ## the non-standard algorithm, that is when type='NS'
  x <- as.matrix(x)  ## makes sure the x is a matrix
  x <- x / sqrt( rowSums(x^2) )  ## makes sure the the data are unit vectors
  n <- nrow(x)  ## sample size
  ina <- as.numeric(ina)
  if ( A >= min(table(ina)) )  A <- min(table(ina)) - 3  ## The maximum 
  ## number  of nearest neighbours to use
  ina <- as.numeric(ina) ## makes sure ina is numeric
  ng <- max(ina)  ## The number of groups
  ## if seed==TRUE then the results will always be the same
  if (seed == TRUE)  set.seed(1234567)
  nu <- sample(1:n, min( n, round(n / M) * M ) )
  ## It may be the case this new nu is not exactly the same
  ## as the one specified by the user
  options(warn = -1)
  mat <- matrix( nu, ncol = M ) # if the length of nu does not fit to a matrix 
  ## a warning message should appear 
  per <- matrix(nrow = M, ncol = A)
  rmat <- nrow(mat)
  dis <- crossprod( t(x) )
  dis <- as.matrix(dis)
  diag(dis) <- 1
  dis <- acos(dis)
  diag(dis) <- 0 
  ## The k-NN algorith is calculated M times. For every repetition a 
  ## fold is chosen and its observations are classified
  for (vim in 1:M) {
    deigma <- mat[, vim]
    id <- ina[ as.vector( t(deigma) ) ]
    ina2 <- ina[ -as.vector( t(deigma) ) ]
    aba <- as.vector( t(deigma) )
    aba <- aba[aba > 0]
    apo <- dis[aba, -aba]
    g <- rep(0, rmat)
    ta <- matrix(nrow = rmat, ncol = ng)
    if (type == "NS") {
      ## Non Standard algorithm
      for (j in 1:A) {
        knn <- j + 1
          for (k in 1:rmat) {
            for (l in 1:ng) {
              dista <- apo[k, ina2 == l]
              if (mesos == TRUE) {
                ta[k, l] <- mean( sort(dista)[1:knn] )
              } else {
                ta[k, l] <- knn / sum( 1/sort(dista)[1:knn] )
              }
            }
          }
        g <- apply(ta, 1, which.min)
        per[vim, j] <- sum(g == id) / rmat
      }
    } else {
      ## Standard algorithm
      for (j in 1:A) {
        knn <- j + 1
        for (k in 1:rmat) {
          xa <- cbind(ina2, apo[k, ])
          qan <- xa[order(xa[, 2]), ]
          sa <- qan[1:knn, 1]
          tab <- table(sa)
          g[k] <- as.integer( names(tab)[ which.max(tab) ] )
        }
        per[vim, j] <- sum(g == id) / rmat
      }
    }
  }
  ela <- colMeans(per)
  bias <- per[ , which.max(ela)] - apply(per, 1, max)  ## TT estimate of bias
  estb <- mean( bias )  ## TT estimate of bias
  names(ela) <- paste("k=", 2:c(A + 1), sep = "")
  plot(2:c(A + 1), ela, type = "b", xlab = "k nearest neighbours", 
  pch = 9, ylab = "Estimated percentage of correct classification")
  percent <- c( max(ela) + estb, estb)
  names(percent) <- c("Estimated percentage", "Estimated bias")
  list( per = ela, percent = percent )  
}