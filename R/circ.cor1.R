################################
#### Cicrular correlation for cicular data I
#### Tsagris Michail 1/2014
#### mtsagris@yahoo.gr
#### References: S Rao Jammalamadaka and A SenGupta (2001)
#### Topics in circular statistics
################################
circ.cor1 <- function(theta, phi, rads = FALSE) {
  ## theta and phi are angular data in degrees or radians
  ## by default they are in degrees
  n <- length(theta)  ## sample size
  ## if the data are in degrees we transform them into radians
  if ( !rads ) {
    theta <- theta * pi/180
    phi <- phi * pi/180
  }
  Rfast2::circ.cor1(theta, phi, pvalue = TRUE)
}
