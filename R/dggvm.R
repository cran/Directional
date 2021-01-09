dggvm <- function(x, param, rads = FALSE, logden = FALSE){

  if (!rads)  x <- x/180 * pi
  n <- length(x)
  z <- param[1]
  k <- param[2]
  m <- param[3]
  a <- param[4]
  phia <- x - a
  ma <- m - a
  cospha <- cos(phia)
  sinpha <- sin(phia)
  cosma <- cos(ma)
  sinma <- sin(ma)
  nzphia <- sqrt(cospha^2 + sinpha^2/z^2)
  nzma <- sqrt(cosma^2 + sinma^2/z^2)
  coszphia <- cospha/nzphia
  sinzphia <- sinpha/(z * nzphia)
  coszma <- cosma/nzma
  sinzma <- sinma/(z * nzma)
  den <-  -k * (coszphia * coszma + sinzphia * sinzma) + 2 * log(nzphia) + log(2 * pi * z) + 
           log(besselI(k, 0, expon.scaled = TRUE)) + k
  
  if ( !logden )  den <- exp(den)
  den
}