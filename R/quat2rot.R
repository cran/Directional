### Code for converting unsigned unit quaternion in R4 to (3 x 3) rotation matrix on SO(3)  ###
quat2rot <-function(x) {
  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  x4 <- x[4]
  X <- matrix(0, 3, 3)
  X[1, 1] <- x1^2 + x4^2 - x2^2 - x3^2
  X[2, 1] <- 2 * ( x3 * x4 + x1 * x2 )
  X[3, 1] <- 2 * ( x1 * x3 - x2 * x4 )
  X[1, 2] <- 2 * ( x1 * x2 - x3 * x4 )
  X[2, 2] <- x2^2 + x4^2 - x1^2 - x3^2
  X[3, 2] <- 2 * ( x2 * x3 + x1 * x4 )
  X[1, 3] <- 2 * ( x2 * x4 + x1 * x3 )
  X[2, 3] <- 2 * ( x2 * x3 - x1 * x4 )
  X[3, 3] <- x3^2 + x4^2 - x1^2 - x2^2
  X
}



