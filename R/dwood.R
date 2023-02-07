################################
#### Wood's bimodal distribution on the sphere
#### Tsagris Michail 1/2016
#### mtsagris@yahoo.gr
#### References: Andrew T.A. Wood (1982), JRSSC, 31(1): 52-58
#### A bimodal distribution on the sphere
################################
dwood <- function(y, param, logden = FALSE) {
  ## y is a two column matrix, where the first column is the latitude and
  ## the second is the longitude, all expressed in degrees
  y <- as.matrix(y)
  if   ( dim(y)[2] == 1 )   y <- t(y)

  y[, 1] <- 90 - y[, 1] ## we want the co-latitude
  y <- y / 180 * pi
  siny1 <- sin( y[, 1] )
  x <- cbind( siny1 * cos(y[, 2]), siny1 * sin(y[, 2]), cos(y[, 1]) )
  #############
  param[1:4] <- param[1:4] / 180 * pi
  gam <- param[1]  ;  del <- param[2]  ;  a <- param[3]
  b <- param[4]    ;  k <- param[5]

  m1 <- c( cos(gam) * cos(del), cos(gam) * sin(del), - sin(gam) )
  m2 <- c( - sin(del), cos(del), 0 )
  m3 <- c( sin(gam) * cos(del), sin(gam) * sin(del), cos(gam) )

  a1 <- as.vector( x %*% m1 )
  a2 <- as.vector( x %*% m2 )
  a3 <- as.vector( x %*% m3 )
  u <- sum(a3)
  down <- sqrt( 1 - a3^2 )
  v <- ( a1^2 - a2^2 ) / down
  w <- 2 * a1 * a2  / down

  den <-  - log(2 * pi) - log( ( exp(k) - exp(-k) ) / k ) +
      k * ( a3 * cos(a) + ( v * cos(b) + w * sin(b) ) * sin(a) )

  if ( !logden )  den <- exp(den)
  den
}

