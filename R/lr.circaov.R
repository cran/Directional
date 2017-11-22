################################
#### ANOVA for cicular data (Likelihood ratio test)
#### Tsagris Michail 1/2015
#### References: Mardia Kanti V. and Jupp Peter E. (2000)
#### Directional statistics, page 136
################################
lr.circaov <- function(u, ina, rads = FALSE) {
  ## u contains all the circular data in radians or degrees
  ## ina is an indicator variable of each sample
  ina <- as.numeric(ina)
  g <- max(ina)  ## how many groups are there
  ## if the data are in degrees we transform them into radians
  if ( !rads )   u <- u * pi/180
  x <- cbind(cos(u), sin(u))
  rsi <- rowsum(x, ina)
  Ri <- sqrt( Rfast::rowsums(rsi^2) )   ## the resultant length of each group
  ni <- tabulate(ina)
  mi <- rsi / ni
  mi <- mi / sqrt( Rfast::rowsums(mi^2) )  ## mean direction of each group
  m <- Rfast::colmeans(x)
  m <- m / sqrt( sum(m^2) )  ## mean direction based on all the data
  m <- matrix(rep(m, g), nrow = g, byrow = TRUE)
  ## Next we estimate the common concentration parameter kappa
  kapa <- circ.summary(u, rads = TRUE, plot = FALSE, fast = TRUE)$kappa
  ## kapa is the estimated concentration parameter based on all the data
  w <- kapa * sum( Ri * Rfast::rowsums((mi - m)^2) )
  pvalue <- pchisq(w, g - 1, lower.tail = FALSE)
  res <- c(w, pvalue, kapa)
  names(res) <- c('test', 'p-value', 'kappa')
  res
}

