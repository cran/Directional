etoa <- function(x) {
  dm <- dim(x)
  n <- dm[1]   ;   d <- dm[2]
  x2 <- x^2
  phi <- matrix(NA, n, d - 1)
  phi[, 1] <- acos( x[, 1] )
  for ( i in 2:(d - 1)  ) phi[, i] <- acos( x[, i] / sqrt( Rfast::rowsums(x2[, i:d]) ) )
  ep <- which(x[, d - 1] < 0 )
  if ( length(ep) > 0 )  phi[ep, d - 1] <-  2 * pi - x[ep, d - 1]
  phi
}
  
    
  

