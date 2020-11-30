## Constructing a rotation matrix on SO(3) from the three euler angles
eul2rot <- function(theta.12, theta.23, theta.13) {
  x.12 <- theta.12
  x.23 <- theta.23
  x.13 <- theta.13
  a1 <- c( cos(x.13) * cos(x.12), cos(x.13) * sin(x.12), sin(x.13) )
  a2 <- c( -cos(x.23) * sin(x.12) - sin(x.23) * sin(x.13) * cos(x.12),
     cos(x.23) * cos(x.12) - sin(x.23) * sin(x.13) * sin(x.12),  sin(x.23) * cos(x.13) )
  a3 <-c( sin(x.23) * sin(x.12) - cos(x.23) * sin(x.13) * cos(x.12),
    -sin(x.23) * cos(x.12) - cos(x.23) * sin(x.13) * sin(x.12), cos(x.23) * cos(x.13) )
  matrix( c(a1, a2, a3), ncol = 3 )
}






