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
  if ( !rads )  u <- u/180 * pi
  x <- cbind( cos(u), sin(u) )
  disa <- tcrossprod(x)
  expa <- exp(disa)
  diag(expa) <-  - Inf  ## we do not want the diagonal elements
  con <- 2 * (n - 1) * pi

  funa <- function(h) {
    A <- expa^h
    f <- rowSums( A )
    sum( log(f) ) - n * log( besselI(h, 0) )
  }

  bar <- optimize(funa, c(low, up), maximum = TRUE)
  res <- c( bar$maximum, bar$objective/n - log(con) )
  names(res) <- c("Optimal h", "cv")
  res
}
