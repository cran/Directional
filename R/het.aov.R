################################
#### ANOVA for hyper-spherical data (Heterogeneous case, kappas not equal)
#### Tsagris Michail 1/2015
#### mtsagris@yahoo.gr
#### References: Mardia Kanti V. and Jupp Peter E. (2000)
#### Directional statistics, page 228
################################
het.aov <- function(x, ina) {
  ## x contains all the data
  ## ina is an indicator variable of each sample
  ina <- as.numeric(ina)
  g <- max(ina)  ## how many groups are there
  p <- dim(x)[2]  ## dimensionality of the data
  ni <- tabulate(ina)  ## group sample sizes
  kapa <- numeric(g)
  mi <- rowsum(x, ina) / ni
  for (i in 1:g)  kapa[i] <- Directional::vmf.mle( x[ina == i, ], fast = TRUE )$kappa
  tw <- Rfast::colsums(kapa * ni * mi)
  Tt <- 2 * ( sum( kapa * ni * sqrt( Rfast::rowsums(mi^2) ) ) - sqrt( sum(tw^2) ) )
  pvalue <- pchisq(Tt, (p - 1) * (g - 1), lower.tail = FALSE)
  res <- c(Tt, pvalue)
  names(res) <- c('test', 'p-value')
  res
}
