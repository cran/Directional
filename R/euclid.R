################################
#### Change between latitude/longitude and Cartesian coordinates
#### Tsagris Michail 1/2014
#### mtsagris@yahoo.gr
################################
euclid <- function(u) {
  u <- as.matrix(u)
  if ( dim(u)[2] == 1 )  u <- t(u)
  u <- pi * u / 180           ## degrees to radians
  clat <- cos(u[, 1])
  U <- cbind(clat * cos(u[, 2]), clat * sin(u[, 2]), sin(u[, 1]))
  colnames(U) <- c("x", "y", "z")
  U
}
