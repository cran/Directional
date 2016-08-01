
################################
#### Kernel density estimation of circular data with a von Mises kernel
#### Tsagris Michail 2/2015
#### mtsagris@yahoo.gr
#### Tuning the bandwidth
################################

vmkde.tune <- function(u, low = 0.1, up = 1, rads = TRUE) {
  ## u is the data
  n <- length(u)  ## sample size

  ## if the data are in degrees we transform them into radians
  if (rads == FALSE)  u <- u/180 * pi
  x <- cbind( cos(u), sin(u) )
  disa <- tcrossprod(x)
  diag(disa) <- 1
  expa <- exp(disa)
  diag(expa) <- NA  ## we do not want the diagonal elements

   funa <- function(h) {
    A <- expa^( 1 / h^2 )
    f <- rowSums( A, na.rm = TRUE )/( (n - 1) * 2 * pi * besselI(1/h^2, 0) ) 
    sum( log(f) ) / n
   }

  bar <- optimize(funa, c(low, up), maximum = TRUE)
  res <- c( bar$maximum, bar$objective )
  names(res) <- c("Optimal h", "cv")
  res
}
