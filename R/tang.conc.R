################################
#### ANOVA for cicular data (Tangential approach for equality of concentration parameters)
#### Tsagris Michail 1/2015
#### mtsagris@yahoo.gr
#### References: Mardia Kanti V. and Jupp Peter E. (2000)
#### Directional statistics, page 141
################################

tang.conc <- function(u, ina, rads = FALSE) {
  ## u contains all the circular data in radians or degrees
  ## ina is an indicator variable of each sample
  n <- length(u)  ## sample size
  ina <- as.numeric(ina)
  ni <- tabulate(ina)
  g <- max(ina)  ## how many groups are there
  ## if the data are in degrees we transform them into radians
  if ( !rads )   u <- u * pi/180
  d <- NULL  ## will store the absolute sinus centred data here
  d2 <- dmi <- numeric(g)
  x1 <- cos(u)
  x2 <- sin(u)
  C <- rowsum(x1, ina)
  S <- rowsum(x2, ina)
  mi <- atan(S/C) + pi * as.numeric(C<0)
  for (i in 1:g) {
    b <- abs( sin( u[ ina == i ] - mi[i] ) )
    d <- c(d, b)
  }

  for (i in 1:g) {
    dmi[i] <- mean( d[ ina == i ] )
    d2[i] <- sum( ( d[ ina == i ] - dmi[i] )^2 )
  }

  mdm <- mean(d)
  Ft <- (n - g) * sum(ni * (dmi - mdm)^2) / ( (g - 1) * sum(d2) )
  pvalue <- pf(Ft, g - 1, n - g, lower.tail = FALSE)
  res <- c(Ft, pvalue)
  names(res) <- c('test', 'p-value')
  res
}
