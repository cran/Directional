################################
#### ANOVA for cicular data (Likelihood ratio test)
#### Tsagris Michail 1/2015
#### References: Mardia Kanti V. and Jupp Peter E. (2000)
#### Directional statistics, page 136
################################

lr.circaov <- function(u, ina, rads = FALSE) {
  ## u contains all the circular data in radians or degrees
  ## ina is an indicator variable of each sample

  n <- length(u)  ## sample size
  ina <- as.numeric(ina)
  g <- max(ina)  ## how many groups are there

  ## if the data are in degrees we transform them into radians
  if ( rads == F )  u <- u * pi/180
  x <- cbind(cos(u), sin(u))
  Ci <- rowsum(x[, 1], ina)
  Si <- rowsum(x[, 2], ina)
  Ri <- sqrt(Ci^2 + Si^2)  ## the resultant length of each group

  ni <- as.vector( table(ina) )
  mi <- rowsum(x, ina) / ni
  mi <- mi/sqrt(rowSums(mi^2))  ## mean direction of each group

  m <- as.vector( Rfast::colmeans(x) )
  m <- m/sqrt(sum(m^2))  ## mean direction based on all the data
  m <- matrix(rep(m, g), nrow = g, byrow = T)

  R <- sqrt( as.vector( sum( Rfast::colsums(x)^2 ) ) )  ## the resultant length based on all the data

  ## Next we estimate the common concentration parameter kappa
  kappa <- circ.summary(u, rads = T, plot = F)$kappa
  ## kappa is the estimated concentration parameter based on all the data

  w <- kappa * sum(Ri * rowSums((mi - m)^2))
  pvalue <- 1 - pchisq(w, g - 1)
  res <- c(w, pvalue, kappa)
  names(res) <- c('test', 'p-value', 'kappa')
  res

}
