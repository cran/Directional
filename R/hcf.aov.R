hcf.aov <- function(x, ina, fc = TRUE) {
  ## x contains all the data
  ## ina is an indicator variable of each sample

  ina <- as.numeric(ina)
  g <- max(ina)  ## how many groups are there
  x <- as.matrix(x)
  x <- x / sqrt( Rfast::rowsums(x^2) )  ## makes sure x are unit vectors
  p <- dim(x)[2]
  n <- dim(x)[1]  ## dimensionality and sample size of the data

  S <- rowsum(x, ina)
  Ri <- sqrt( Rfast::rowsums(S^2) )  ## the resultant length of each group
  S <- Rfast::colsums(x) 
  R <- sqrt( sum(S^2) )  ## the resultant length based on all the data

  ## Next we stimate the common concentration parameter kappa
  kappa <- vmf(x)$kappa
  ## kappa is the estimated concentration parameter based on all the data

  Ft <- (n - g) * (p - 1) * (sum(Ri) - R) /( (g - 1) * (p - 1) * (n - sum(Ri)) )

  if (fc == TRUE) {  ## correction is used
    if (p == 3) {
      Ft <- kappa * (1/kappa - 1/(5 * kappa^3)) * Ft
    } else if (p > 3)  {
      Ft <- kappa * ( 1/kappa - (p - 3)/(4 * kappa^2) - (p - 3)/(4 * kappa^3) ) * Ft
    }
  }

  pvalue <- pf(Ft, (g - 1) * (p - 1), (n - g) * (p - 1), lower.tail = FALSE)
  res <- c(Ft, pvalue, kappa)
  names(res) <- c('test', 'p-value', 'kappa')
  res

}
