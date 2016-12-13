###############################
#### Linear cicrular correlation
#### Tsagris Michail 3/2014
#### mtsagris@yahoo.gr
#### References: Mardia Kanti V. and Jupp Peter E. (2000)
#### Directional statistics
################################

circlin.cor <- function(theta, x, rads = FALSE) {
  ## theta is a angular variable in degrees by default
  ## x is euclidean variable or a matrix containing euclidean variables

  n <- length(theta)  ## sample size
  if ( !rads )  theta <- theta * pi/180
  costheta <- cos(theta)
  sintheta <- sin(theta)

  rxc <- as.numeric( cor(costheta, x) )  ## and cos(theta) and x correlation
  rxs <- as.numeric( cor(sintheta, x) )  ## sin(theta) and x correlation
  rcs <- cor(costheta, sintheta)  ## cos(theta) and sin(theta) correlation
  R2xt <- (rxc^2 + rxs^2 - 2 * rxc * rxs * rcs)/(1 - rcs^2)

  ## linear-circular correlation
  Ft <- (n - 3) * R2xt / (1 - R2xt)  ## F-test statistic value
  pvalue <- pf(Ft, 2, n - 3, lower.tail = FALSE)
  res <- cbind(R2xt, pvalue)
  colnames(res) <- c('R-squared', 'p-value')
  res

}
