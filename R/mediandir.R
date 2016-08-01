################################
#### Median direction
#### Tsagris Michail 1/2016
#### mtsagris@yahoo.gr
#### References: Fisher, N. I. (1985). Spherical medians.
#### Journal of the Royal Statistical Society. Series B, 47(2): 342-348.
#### Fisher, N. I., Lewis, T., & Embleton, B. J. (1987).
#### Statistical analysis of spherical data. Cambridge university press.
#### Cabrera, J., & Watson, G. S. (1990). On a spherical median related distribution.
#### Communications in Statistics-Theory and Methods, 19(6): 1973-1986.
################################

mediandir = function(x) {
  ## x is the directional data
  x <- as.matrix(x)
  x <- x / sqrt( rowSums(x^2) )
  n <- nrow(x)
  p <- ncol(x)

  pa <- as.vector( Rfast::colMedians(x) )
  u1 <- pa / sqrt( sum(pa^2) )
  ww <- as.vector( sqrt( 1 - ( x %*% u1 )^2 ) )
  u2 <- as.vector( Rfast::colsums(x / ww ) )
  u2 <- u2 / sqrt( sum( u2^2 ) )

  i <- 2
  while ( sum( abs (u2 - u1 ) ) > 1e-10 ) {
    i <- i +1
    u1<- u2
    ww <- as.vector( sqrt( 1 - ( x %*% u1 )^2 ) )
    u2 <- as.vector( Rfast::colsums (x / ww ) )
    u2 <- u2 / sqrt( sum( u2^2 ) )
  }

  u2
}
