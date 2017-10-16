###############################
#### Linear cicrular correlation
#### Tsagris Michail 3/2014
#### mtsagris@yahoo.gr
#### References: Mardia Kanti V. and Jupp Peter E. (2000)
#### Directional statistics
################################
circlin.cor <- function(theta, x, rads = FALSE) {
  ## theta is a angular variable in degrees by default
  ## x is euclidean variable or a matrix containing euclidean variables
  if ( !rads )  theta <- theta * pi/180
  Rfast::circlin.cor(theta, x)
}
