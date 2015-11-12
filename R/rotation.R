################################
#### Rotation matrix 
#### Tsagris Michail 10/2013 
#### mtsagris@yahoo.gr
#### References: G. J. A. Amaral, I. L. Dryden & Andrew T. A. Wood (2007, JASA)
#### Pivotal Bootstrap Methods for k-Sample Problems in Directional Statistics and Shape Analysis. 
################################

rotation <- function(a, b) {
  ## a and b are two unit vectors
  ## Calculates the rotation matrix
  ## to move a to b along the geodesic path
  ## on the unit sphere which connects a to b
  p <- length(a)
  c <- a - b * (a %*% t( t(b)) )
  c <- c/sqrt( sum(c^2) )
  A <- t( t(b) ) %*% c - t( t(c) ) %*% b
  theta <- acos(sum(a * b))
  diag(p) + sin(theta) * A + (cos(theta) - 1) * ( t( t(b) ) %*% b + t( t(c) ) %*% c )
 }