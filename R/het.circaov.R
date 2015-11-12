################################
#### ANOVA for cicular data (Heterogeneous case, kappas not equal)
#### Tsagris Michail 1/2015
#### mtsagris@yahoo.gr
#### References: Mardia Kanti V. and Jupp Peter E. (2000)
#### Directional statistics, page 141-142
################################

het.circaov <- function(u, ina, rads = F) {
  ## u contains all the circular data in radians or degrees
  ## ina is an indicator variable of each sample
  n <- length(u)  ## sample size
  ina <- as.numeric(ina)
  g <- max(ina)  ## how many groups are there
  ni <- as.vector(table(ina))
  ## if the data are in degrees we transform them into radians
  if (rads == F) u <- u * pi/180 
  kappa <- numeric(g)
  x1 <- cos(u)
  x2 <- sin(u)
  C <- aggregate(x1, by = list(ina), sum)
  S <- aggregate(x2, by = list(ina), sum)
  C <- as.vector( C[, -1] )
  S <- as.vector( S[, -1] )
  mi <- atan(S/C) + pi * as.numeric(C<0)
  Ri <- sqrt(C^2 + S^2)  ## the resultant length of each group
  for (i in 1:g) kappa[i] <- circ.summary(u[ina == i], rads = T, plot = F)$kappa
  ## kappa contains the estimated concentration parameters of each group
  Rw <- ( (sum(kappa * Ri * cos(mi)) )^2 + ( sum(kappa * Ri * sin(mi)))^2 )^0.5
  T <- 2 * (sum(kappa * Ri) - Rw)
  pvalue <- 1 - pchisq(T, g - 1)
  res <- c(T, pvalue)
  names(res) <- c('test', 'p-value')
  res
}