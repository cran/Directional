################################
#### Kernel density estimation of directional data with a von Mises kernel
#### Tuning the bandwidth
#### Tsagris Michail 8/2015
#### mtsagris@yahoo.gr
#### Tuning the bandwidth
################################

vmfkde.tune <- function(x, low = 0.1, up = 1) {
  ## x is the data
  p <- dim(x)[2]  ## dimensionality of the data
  n <- dim(x)[1]  ## sample size of the data
  d <- tcrossprod( x )
  diag(d) <- NA  ##  we do not want to take the diagonal elements
  con <- (2 * pi)^(p/2)

   funa <- function(h) {
    A <- d/h^2
    cpk <- (1/h^2)^(p/2 - 1) / con / besselI(1/h^2, p/2 - 1)
    f <- rowSums( exp(A + log(cpk)), na.rm = TRUE )/(n - 1)
    mean( log(f) )
  }

  a <- optimize(funa, c(low, up), maximum = TRUE)
  res <- c(a$maximum, a$objective)
  names(res) <- c("Optimal h", "cv")
  res
}
