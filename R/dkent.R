dkent <- function(y, G, param, logden = FALSE ) {
  k <- param[1]
  b <- param[2]
  gam <- c(0, k, 0)
  lam <- c(0, -b, b)
  ckb <- Directional::fb.saddle(gam, lam)[3]
  y <- as.matrix(y)
  if ( dim(y)[2] == 1 )  y <- t(y)
  den <-  -ckb + k * y %*% G[, 1] + b * (y %*% G[, 2])^2 - b * (y %*% G[, 3])^2
  den <- as.vector(den)
  if ( !logden )  den <- exp(den)
  den
}
