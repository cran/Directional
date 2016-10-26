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
  if ( rads == FALSE ) {
    theta <- theta * pi/180
    phi <- phi * pi/180
  }

  ## We calculate the mean of each vector
  m1 <- circ.summary( theta, rads = TRUE, plot = FALSE )$mesos
  m2 <- circ.summary( phi, rads = TRUE, plot = FALSE )$mesos
  sintheta <- sin(theta - m1)
  sinphi <- sin(phi - m2)

  up <- sum( sintheta * sinphi )
  down <- sqrt( sum( sintheta ^2 ) * sum( sinphi^2 ) )
  rho <- up/down  ## circular correlation
  lam22 <- sum( sintheta^2 * sinphi^2 ) / n
  lam02 <- sum( sinphi^2 ) / n
  lam20 <- sum( sintheta^2 ) / n
  zrho <- sqrt(n) * sqrt( lam02 * lam20/lam22 ) * rho
  pvalue <- 2 * pnorm( abs(zrho), lower.tail = FALSE )

  res <- c(rho, pvalue)
  names(res) <- c("rho", "p-value")
  res

}
