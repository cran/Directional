################################
#### Change between latitude/longitude and Cartesian coordinates
#### Tsagris Michail 1/2014
#### mtsagris@yahoo.gr
################################
euclid.inv <- function(U) {
  U <- as.matrix(U)
  if ( dim(U)[2] == 1 )  U <- t(U)
  lat <- asin(U[, 3])
  lon <- atan2(U[, 2], U[, 1])
  res <- 180 * cbind(lat, lon) / pi   ## back to degrees
  colnames(res) <- c("latitude", "longitude")
  res
}
