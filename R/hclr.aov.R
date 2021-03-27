hclr.aov <- function(x, ina) {

  ina <- as.numeric(ina)
  g <- max(ina)  ## how many groups are there
  p <- dim(x)[2]  ## dimensionality of the data
  n <- dim(x)[1]  ## sample size of the data
  S <- rowsum(x, ina)
  Ri <- sqrt( Rfast::rowsums(S^2) )  ## the resultant length of each group
  S <- Rfast::colsums(S)
  R <- sqrt( sum(S^2) )  ## the resultant length based on all the data

  Apk <- function(p, k)  besselI(k, p/2, expon.scaled = TRUE) / besselI(k, p/2 - 1, expon.scaled = TRUE)

  ## Next we estimate the common concentration parameter kappa under H0 and H1
  Rk <- R/n
  k1 <- Rk * (p - Rk^2)/(1 - Rk^2)
  k2 <- k1 - (Apk(p, k1) - Rk) / ( 1 - Apk(p, k1)^2 - (p - 1)/k1 * Apk(p, k1) )
  while ( abs(k2 - k1) > 1e-07 ) {
    k1 <- k2
    k2 <- k1 - (Apk(p, k1) - Rk) / (1 - Apk(p, k1)^2 - (p - 1)/k1 * Apk(p, k1) )
  }
  k0 <- k2  ## concentration parameter under H0

  Rk <- sum(Ri)/n
  k1 <- Rk * (p - Rk^2)/(1 - Rk^2)
  k2 <- k1 - (Apk(p, k1) - Rk) / ( 1 - Apk(p, k1)^2 - (p - 1)/k1 * Apk(p, k1) )

  while ( abs(k2 - k1) > 1e-07 ) {
    k1 <- k2
    k2 <- k1 - (Apk(p, k1) - Rk) / ( 1 - Apk(p, k1)^2 - (p - 1)/k1 * Apk(p, k1) )
  }
  k1 <- k2  ## concentration parameter under H1

  I0 <- besselI( k0, p/2 - 1 )
  I1 <- besselI( k1, p/2 - 1 )
  Apk0 <- Apk(p, k0)
  Apk1 <- Apk(p, k1)
  up <- k0^(p/2 - 1) / I0 * exp(k0 * Apk0)
  down <- k1^(p/2 - 1) / I1 * exp(k1 * Apk1)
  P <- (n - g) / (g - 1) * ( (up / down)^( -2/(p - 1) ) - 1 )
  pvalue <- pf(P, (g - 1) * (p - 1), (n - g) * (p - 1), lower.tail = FALSE)
  res <- c(P, pvalue)
  names(res) <- c('P', 'p-value')
  res
}

