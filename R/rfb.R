################################
#### Simulating from a Fisher-Bingham distribution
#### Tsagris Michail 03/2019
#### mtsagris@yahoo.gr
#### References: A new method to simulate the Bingham and related distributions in
#### directional data analysis with applications
#### Kent J.T., Ganeiber A.M. and Mardia K.V. (2013)
#### http://arxiv.org/pdf/1310.8110v1.pdf
################################

rfb <- function(n, k, m, A) {
  ## n is the required sample size
  ## k is the concentration parameter, the Fisher part
  ## m is the mean vector, the Fisher part
  ## A is the symmetric matrix, the Bingham part
  m <- m / sqrt( sum(m^2) )
  m0 <- c(0, 1, 0)
  mu <- c(0, 0, 0)
  B <- rotation(m0, m)
  q <- length(m0)
  A1 <- A + k/2 * ( diag(q) - m0 %*% t(m0) )
  eig <- eigen(A1)
  lam <- eig$values
  V <- eig$vectors
  lam <- lam - lam[q]
  lam <- lam[-q]
  x <- f.rbing(n, lam, fast = TRUE)$X  ## Chris and Theo's code
  x <- tcrossprod(x, V) ## simulated data
  u <- log( Rfast2::Runif(n) )
  ffb <- k * x[, 2]  - Rfast::rowsums( x %*% A * x )
  fb <- k - Rfast::rowsums( x %*% A1 * x )
  x1 <- x[u <= c(ffb - fb), ]
  n1 <- dim(x1)[1]

  while (n1 < n) {
    x <- Directional::f.rbing(n - n1, lam, fast = TRUE)$X  ## Chris and Theo's code
    x <- tcrossprod(x, V) ## simulated data
    u <- log( runif(n - n1) )
    ffb <- k * x[, 2]  - Rfast::rowsums( x %*% A * x )
    fb <- k - Rfast::rowsums( x %*% A1 * x )
    x1 <- rbind(x1, x[u <= c(ffb - fb), ])
    n1 <- dim(x1)[1]
  }

  tcrossprod(x1, B) ## simulated data with the wanted mean direction
}



################################
#### Simulating from a Fisher-Bingham distribution
#### Tsagris Michail 05/2014
#### mtsagris@yahoo.gr
#### References: A new method to simulate the Bingham and related distributions in
#### directional data analysis with applications
#### Kent J.T., Ganeiber A.M. and Mardia K.V. (2013)
#### http://arxiv.org/pdf/1310.8110v1.pdf
################################

# rfb <- function(n, k, m, A) {
#   ## n is the required sample size
#   ## k is the concentration parameter, the Fisher part
#   ## m is the mean vector, the Fisher part
#   ## A is the symmetric matrix, the Bingham part
#   m <- m / sqrt( sum(m^2) )
#   m0 <- c(0, 1, 0)
#   mu <- c(0, 0, 0)
#   B <- rotation(m0, m)
#   q <- length(m0)
#   A1 <- A + k/2 * ( diag(q) - m0 %*% t(m0) )
#   eig <- eigen(A1)
#   lam <- eig$values
#   V <- eig$vectors
#   lam <- lam - lam[q]
#   lam <- lam[-q]
#   x1 <- matrix( 0, n, 3 )
#   i <- 1
#
#   while (i <= n) {
#     x <- f.rbing(1, lam, fast = TRUE)$X  ## Chris and Theo's code
#     x <- tcrossprod(x, V) ## simulated data
#     u <- log( runif(1) )
#     ffb <- k * x[, 2]  - sum( x %*% A * x )
#     fb <- k - sum( x %*% A1 * x )
#     if ( u <= c(ffb - fb) ) {
#       x1[i, ] <- x
#       i <- i + 1
#     }
#   }
#
#   tcrossprod(x1, B) ## simulated data with the wanted mean direction
# }
