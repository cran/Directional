## Probability density function of the
## von Mises-Fisher distribution
## May 2016
## References: Arthur Pewsey, Markus Neuhauser, and Graeme D. Ruxton (2013)
## Circular Statistics in R


pvm <- function(theta, m, k, rads = FALSE) {

  if ( rads == FALSE )  {
     u <- u * pi / 180
     m <- m * pi / 180
  }

   theta <- theta %% (2 * pi)

   if ( k > 0 ) {
     theta <- theta %% (2 * pi)
     f <- 2 * pi * besselI(k, 0)
     funa <- function(u)  exp(k * cos(u - m) )
     prob <- as.numeric( integrate(funa, 0, theta)$value ) / f

   } else {
     prob <- theta / ( 2 * pi )
   }

  prob

}
