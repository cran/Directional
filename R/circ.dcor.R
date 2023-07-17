circ.dcor <- function(theta, phi, rads = FALSE) {
  if ( !rads ) {
    theta <- theta * pi/180
    phi <- phi * pi/180
  }
    
  Rfast::dcor( cbind( cos(theta), sin(theta) ), cbind( cos(phi), sin(phi) ) )
}
