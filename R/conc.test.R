################################
#### ANOVA for cicular data (Test for equality of concentration parameters)
#### Tsagris Michail 2/2015
#### mtsagris@yahoo.gr
#### References: Mardia Kanti V. and Jupp Peter E. (2000)
#### Directional statistics, page 139-140
################################
conc.test <- function(u, ina, rads = FALSE) {
  ## u contains all the circular data in rads or degrees
  ## ina is an indicator variable of each sample
  n <- length(u)  ## sample size
  ina <- as.numeric(ina)
  g <- max(ina)  ## how many groups are there
  ni <- tabulate(ina)
  if ( !rads )  u <- u * pi/180
  ## if the data are in degrees we transform them into radians
  x1 <- cos(u)
  x2 <- sin(u)
  Ci <- Rfast::group(x1, ina, method = "mean")  ## rowsum(x1, ina) / ni
  Si <- Rfast::group(x2, ina, method = "mean")  ## rowsum(x2, ina) / ni
  Rbi <- sqrt( Ci^2 + Si^2 )
  ## Ri is the mean resultant length of each group
  C <- sum( Ci * ni ) / n   ### sum(x1)/n
  S <- sum( Si * ni ) / n   ### sum(x2)/n
  Rb <- sqrt(C^2 + S^2)  ## the mean resultant length of all the data

  if (Rb < 0.45) {
    ## case 1
    g1 <- wi <- numeric(g)
    wi <- 4 * (ni - 4)/3
    g1 <- asin( sqrt(3/8) * 2 * Rbi )
    U <- sum(wi * g1^2) - ( sum(wi * g1) )^2 / sum(wi)
    pvalue <- pchisq(U, g - 1, lower.tail = FALSE)
    mess <- paste('The mean resultant length is less than 0.45. U1 was calculated')

  } else if (Rb >= 0.45 & Rb <= 0.7) {
    ## case 2
    g2 <- wi <- numeric(g)
    wi <- (ni - 3) / 0.798
    g2 <- asinh( (Rbi - 1.089) / 0.258 )
    U <- sum( wi * g2^2) - (sum(wi * g2) )^2 / sum(wi)
    pvalue <- pchisq(U, g - 1, lower.tail = FALSE)
    mess <- paste('The mean resultant length is between 0.45 and 0.7. U2 was calculated')

  } else  if (Rb > 0.7) {
    ## case 3
    Ri <- Rbi * ni
    vi <- ni - 1
    v <- n - g
    d <-  ( sum(1/vi) - 1/v ) / ( 3 * (g - 1) )
    U <- ( v * log( (n - sum(Ri))/v ) - sum( vi * log((ni - Ri)/vi) ) ) / (1 + d)
    pvalue <- pchisq(U, g - 1, lower.tail = FALSE)
    mess <- paste('The mean resultant length is more than 0.7. U3 was calculated')
  }

  res <- c(U, pvalue)
  names(res) <- c('test', 'p-value')
  list(mess = mess, res = res)
}
