################################
#### Simulation of random rotation matrices
#### Matrices in SO(p)
################################

rsop <- function(n, p) {

  a <- c(1, numeric(p - 1) )
  A <- array( dim = c(p, p, n) )
  Ip <- diag(p)

  for (i in 1:n) {
    b <- rnorm(p)
    b <- b / sqrt( sum(b^2) )
    ca <- a - b * b[1]
    ca <- ca / sqrt( sum(ca^2) )
    B <- tcrossprod(b, ca)
    B <- B - t(B)
    theta <- acos( b[1] )
    A[, , i] <- Ip + sin(theta) * B + ( cos(theta) - 1 ) * ( tcrossprod(b) + tcrossprod(ca) )
  }

  if (n == 1)  A <- as.matrix( A[, , 1] )
  A
}
