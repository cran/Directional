################################
#### Kernel density estimation of circular data with a von Mises kernel
#### Tsagris Michail 2/2015 
#### mtsagris@yahoo.gr
#### Tuning the bandwidth
################################

vmkde.tune <- function(u, h = seq(0.1, 1, by = 0.01), rads = T, 
  plot = TRUE, ncores = 4) {
  ## u is the data
  ## h is the bandwidth grid of values
  ## nc is the number of cores you want to use
  ## requires(doParallel)
  requireNamespace("doParallel", quietly = TRUE)  
  n <- length(u)  ## sample size
  ## if the data are in degrees we transform them into radians
  if (rads == FALSE) u <- u/180 * pi  
  disa <- dist(u, diag = T, upper = T)
  disa <- as.matrix(disa)
  options(warn = -1)
  nc <- ncores
  val <- matrix(h, ncol = nc) ## if the length of h is not equal to the 
  ## dimensions of the matrix val a warning message should appear 
  ## but with options(warn = -1) you will not see it
  cl <- makePSOCKcluster(nc)
  registerDoParallel(cl)
  j = NULL
  ww <- foreach(j = 1:nc, .combine = cbind) %dopar% {
   ba <- val[, j]
    for (l in 1:length(val[, j])) {
      A <- exp( cos(disa) / ba[l]^2 )
      diag(A) <- NA  ## we do not want the diagonal elements
      f <- rowSums(A, na.rm = T)/((n - 1) * 2 * pi * besselI( 1/ba[l]^2, 0) )
      ba[l] <- mean( log(f) )
    }
   return(ba)
  }
  stopCluster(cl)
  cv <- as.vector(ww)[1:length(h)] 
  if (plot == TRUE) {
    plot(h, cv, type = "l")
  }
  list(hopt = h[which.max(cv)], cv = cv)
}