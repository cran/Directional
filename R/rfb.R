################################
#### Simulating from a Fisher-Bingham distribution 
#### Tsagris Michail 05/2014 
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
  m <- m / sqrt(sum(m^2))
  m0 <- c(0, 1, 0)
  B <- rotation(m0, m)
  q <- length(m0)
  A1 <- A + k/2 * (diag(q) - m0 %*% t(m0))
  eig <- eigen(A1)
  lam <- eig$values
  V <- eig$vectors
  lam <- lam - lam[q]
  lam <- lam[-q]
  x <- f.rbing(5 * n, lam)$X  ## Chris and Theo's code
  x <- x %*% t(V)
  u <- log(runif(5 * n))
  ffb <- k * x %*% m0 - diag(x %*% A %*% t(x))
  fb <- k - diag(x %*% A1 %*% t(x))
  x1 <- x[u <= c(ffb - fb), ]
  x <- x1[sample(1:nrow(x1), n), ]
  x %*% t(B)  ## simulated data with the wanted mean direction
}