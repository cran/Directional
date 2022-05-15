################################
#### Discrminant analysis for directional data
#### assuming a von Mises-Fisher distribution
#### Cross-validation for the performance
#### Tsagris Michail 03/2014
#### mtsagris@yahoo.gr
#### References: J. E. Morris and P. J. Laycock (1974)
#### Discriminant Analysis of Directional Data (Biometrika)
################################
vmf.da <- function(x, ina, fraction = 0.2, R = 200, seed = NULL) {
  ## x is the data set
  ## ina is the group indicator variable
  ## fraction denotes the percentage of the sample to be used as the test sample
  ## R is the number of cross validations
  runtime <- proc.time()
  p <- dim(x)[2]  ## p is the dimensionality of the data
  per <- numeric(R)
  n <- dim(x)[1]  ## sample size
  ina <- as.numeric(ina)
  frac <- round(fraction * n)
  g <- max(ina)  ## how many groups are there
  mesi <- matrix(nrow = g, ncol = p)
  k <- numeric(g)
  ## if seed==TRUE then the results will always be the same
  if ( !is.null(seed) )  set.seed(seed)

  for (i in 1:R) {
    nu <- Rfast2::Sample.int(n, frac)
    mod <- Rfast::multivmf.mle(x[-nu, ], ina[-nu], ell = FALSE)
    ki <- mod$ki
    mat <- (p/2 - 1) * log(ki) + ki * tcrossprod(mod$mi, x[nu, ]) - log( besselI(ki, p/2 - 1, expon.scaled = TRUE) ) - ki
    est <- Rfast::colMaxs(mat)
    per[i] <- sum( est == ina[nu] ) / frac
  }

  percent <- mean(per)
  s1 <- sd(per)
  s2 <- sqrt( percent * (1 - percent) / R )
  conf1 <- c(percent - 1.96 * s1, percent + 1.96 * s1)  ## 1st way of a CI
  conf2 <- c(percent - 1.96 * s2, percent + 1.96 * s2)  ## 2nd way of a CI
  ## next we check if the confidence limits exceeds the allowed limits
  if (conf1[2] > 1) conf1[2] <- 1
  if (conf1[1] < 0) conf1[1] <- 0
  if (conf2[2] > 1) conf2[2] <- 1
  if (conf2[1] < 0) conf2[1] <- 0

  conf3 <- quantile(per, probs = c(0.025, 0.975))  ## 3rd way of a CI
  ci <- rbind(conf1, conf2, conf3)
  runtime <- proc.time() - runtime
  colnames(ci) <- c("2.5%", "97.5%")
  rownames(ci) <- c("standard", "binomial", "empirical")
  percent <- c(percent, s1, s2)
  names(percent) <- c('percent', 'sd1', 'sd2')
  list(percent = percent, ci = ci, runtime = runtime)
}
