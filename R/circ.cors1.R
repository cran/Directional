circ.cors1 <- function(theta, phi, rads = FALSE) {
  if ( !rads ) {
    theta <- theta * pi/180
    phi <- phi * pi/180
  }
  Rfast2::circ.cors1(theta, phi, pvalue = TRUE)
}
