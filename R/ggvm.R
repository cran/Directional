################################
#### Geometrically generalised von Mises distribution
#### mtsagris@yahoo.gr
#### References: Fisher, N. I. (1985). Spherical medians.
#### Dietrich, T., & Richter, W. D. (2016).
#### Classes of geometrically generalized von Mises distributions.
#### Sankhya B, 1-39.
################################

ggvm <- function(phi, rads = FALSE) {

  if ( !rads )  phi <- phi / 180 * pi
  n <- length(phi)

  likel <- function(pa, phi) {
    z <- abs( pa[1] )    ;    k <- exp( pa[2] )
    m <- pa[3]    ;    a <- pa[4]

    phia <- phi - a
    ma <- m - a
    cospha <- cos(phia)
    sinpha <- sin(phia)
    cosma <- cos(ma)
    sinma <- sin(ma)

    nzphia <- sqrt( cospha^2 + sinpha^2 / z^2 )
    nzma <- sqrt( cosma^2 + sinma^2 / z^2 )
    coszphia <- cospha / nzphia
    sinzphia <- sinpha / ( z * nzphia )
    coszma <- cosma / nzma
    sinzma <- sinma / ( z * nzma )

    ell <-  - k * sum( coszphia * coszma + sinzphia * sinzma) +
            2 * sum( log(nzphia) ) + n * log(2 * pi * z) +
            n * ( log( besselI(k, 0, expon.scaled = TRUE) ) + k )
    ell
  }

  qa <- optim( c(rnorm(3, 0, 0.1), runif(1, 0, pi) ), likel, phi = phi, control = list(maxit = 5000) )
  qa <- optim( qa$par, likel, phi = phi )
  qa <- optim( qa$par, likel, phi = phi )

  param <- c( abs( qa$par[1]) , exp(qa$par[2]), qa$par[3], qa$par[4] )
  names(param) <- c("Zeta", "kappa", "mu", "alpha")

  list(loglik = -qa$value, param = param)

}




