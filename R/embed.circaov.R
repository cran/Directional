################################
#### ANOVA for cicular data (Embedding approach)
#### Tsagris Michail 1/2015
#### mtsagris@yahoo.gr
#### References: Mardia Kanti V. and Jupp Peter E. (2000)
#### Directional statistics, page 138-139
################################

embed.circaov <- function(u, ina, rads = FALSE) {
  ## u contains all the circular data in radians or degrees
  ## ina is an indicator variable of each sample

  n <- length(u)  ## sample size
  ina <- as.numeric(ina)
  ni <- as.vector(table(ina))
  g <- max(ina)  ## how many groups are there
  ## if the data are in degrees we transform them into radians

  if ( rads == FALSE )  u <- u * pi/180
  x1 <- cos(u)
  x2 <- sin(u)
  Ci <- rowsum(x1, ina)
  Si <- rowsum(x2, ina)
  Rbi <- sqrt( Ci^2 + Si^2 )/ni

  C <- sum(Ci)
  S <- sum(Si)
  Rbar <- sqrt(C^2 + S^2)/n  ## the mean resultant length based on all the data
  ## Next we estimate the common concentration parameter kappa

  kappa <- circ.summary(u, rads = TRUE, plot = FALSE)$kappa
  ## kappa is the estimated concentration parameter based on all the data

  Fb <- ( (sum(ni * Rbi^2) - n * Rbar^2 )/(g - 1) ) /
    ( (n - sum(ni * Rbi^2) )/(n - g) )
  Fc <- ( 1 - 1/(5 * kappa) - 1/(10 * kappa^2) ) * Fb
  pvalue <- pf(Fc, g - 1, n - g, lower.tail = FALSE)
  res <- c(Fc, pvalue, kappa)
  names(res) <- c('test', 'p-value', 'kappa')
  res

}
