lr.perm <- function(x1, x2, B = 999) {

  n1 <- dim(x1)[1]    ;  n2 <- dim(x2)[1]
  x <- rbind(x1, x2)
  ina <- c( rep(1, n1), rep(2, n2) )
  ni <- c(n1, n2)
  p <- dim(x)[2]
  n <- n1 + n2

  S <- rowsum(x, ina)
  Ri <- sqrt( Rfast::rowsums(S^2) )  ## the resultant length of each group
  S <- Rfast::colsums(S)
  R <- sqrt( sum(S^2) )  ## the resultant length based on all the data

  Apk <- function(p, k)  besselI(k, p/2, expon.scaled = TRUE) / besselI(k, p/2 - 1, expon.scaled = TRUE)

  ## Next we stimate the common concentration parameter kappa under H0 and H1
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

  apk0 <- (1 - p/2) * log(k0/2) + lgamma(p/2) + log( besselI( k0, p/2 - 1, expon.scaled = TRUE ) ) + k0
  apk1 <- (1 - p/2) * log(k1/2) + lgamma(p/2) + log( besselI( k1, p/2 - 1, expon.scaled = TRUE ) ) + k1
  w <- k1 * sum(Ri) - k0 * R - n * apk1 + n * apk0

  wp <- numeric(B)
  for (i in 1:B) {
    ind <- sample(ina, n)
    S <- rowsum(x, ind)
    Ri <- sqrt( Rfast::rowsums(S^2) )  ## the resultant length of each group

    Apk <- function(p, k)  besselI(k, p/2, expon.scaled = TRUE) / besselI(k, p/2 - 1, expon.scaled = TRUE)

    ## Next we stimate the common concentration parameter kappa under H0 and H1
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

    apk0 <- (1 - p/2) * log(k0/2) + lgamma(p/2) + log( besselI( k0, p/2 - 1, expon.scaled = TRUE ) ) + k0
    apk1 <- (1 - p/2) * log(k1/2) + lgamma(p/2) + log( besselI( k1, p/2 - 1, expon.scaled = TRUE ) ) + k1
    wp[i] <- k1 * sum(Ri) - k0 * R - n * apk1 + n * apk0
  }

  pvalue <- ( sum(wp > w) + 1 ) / (B + 1)
  res <- c(2 * w, pvalue)
  names(res) <- c('w', 'p-value')
  res
}

