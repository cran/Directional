################################
#### Kernel density estimation of circular data with a von Mises kernel
#### Tsagris Michail 2/2015 
#### mtsagris@yahoo.gr
#### Tuning the bandwidth
################################

vmfkde.tune <- function(x, h = seq(0.1, 1, by = 0.01), plot = TRUE) {
  ## x is the data
  ## h is the bandwidth grid you want
  x <- as.matrix(x)  ## makes sure x is a matrix
  x <- x/sqrt(rowSums(x^2))  ## makes sure x is directional data
  p <- ncol(x)  ## dimensionality of the data
  n <- nrow(x)  ## sample size of the data
  cv <- numeric(length(h))
  d <- crossprod(t(x))
  diag(d) <- NA  ##  we do not want to take the diagonal elements
  for (j in 1:length(h)) {
    A <- d/h[j]^2 
    cpk <- ((1/h[j]^2)^(p/2 - 1))/((2 * pi)^(p/2) * besselI(1/h[j]^2, p/2 - 1))
    f <- rowSums( exp(A + log(cpk)), na.rm = TRUE )/(n - 1)
    cv[j] <- mean(log(f))
  }
  if (plot == TRUE)  plot(h, cv, type = "l")
  list(hopt = h[which.max(cv)], cv = cv)
}