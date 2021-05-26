################################
#### ANOVA for hyper-spherical data (Embedding approach)
#### Tsagris Michail 1/2015
#### mtsagris@yahoo.gr
#### References: Mardia Kanti V. and Jupp Peter E. (2000)
#### Directional statistics, page 225-226
################################
embed.aov <- function(x, ina) {
  ## x contains all the data
  ## ina is an indicator variable of each sample
  ina <- as.numeric(ina)
  g <- max(ina)  ## how many groups are there
  ni <- tabulate(ina)
  p <- dim(x)[2]  ## dimensionality of the data
  n <- dim(x)[1]  ## sample size of the data
  S <- rowsum(x, ina) / ni
  Rbi <- sqrt( Rfast::rowsums(S^2) ) ## the mean resultant length of each group
  S <- Rfast::colmeans(x)
  Rbar <- sqrt( sum(S^2) )  ## the mean resultant length based on all the data
  Ft <- (n - g) * ( sum(ni * Rbi^2) - n * Rbar^2) / ( (g - 1) * ( n - sum(ni * Rbi^2) ) )
  pvalue <- pf(Ft, (g - 1) * (p - 1), (n - g) * (p - 1), lower.tail = FALSE)
  res <- c(Ft, pvalue)
  names(res) <- c('test', 'p-value')
  res
}

