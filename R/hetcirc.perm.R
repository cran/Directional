hetcirc.perm <- function(u1, u2, rads = TRUE, B = 999) {

  u <- c(u1, u2)
  if ( !rads ) u <- u * pi/180
  ina <- c( rep(1, length(u1) ), rep(2, length(u2) ) )
  N <- length(u)
  ni <- tabulate(ina)
  kapa <- numeric(2)
  x1 <- cos(u)
  x2 <- sin(u)
  C <- Rfast::group(x1, ina)
  S <- Rfast::group(x2, ina)
  mi <- atan(S/C) + pi * as.numeric(C < 0)
  Ri <- sqrt(C^2 + S^2)
  ki <- (1.28 - 0.53 * Ri^2/ni^2) * tan(0.5 * pi * Ri/ni)
  coni <- ki

  for (i in 1:2) {
    n <- ni[i]
    coni[i] <- sum( cos( u[ina == i] - mi[i] ) )
    con <- coni[i]
    k1 <- ki[i]
    if (k1 < 710) {
      der <- con - n * besselI(k1, 1, expon.scaled = TRUE)/besselI(k1, 0, expon.scaled = TRUE)
      a <- besselI(k1, 0)^2/2 + besselI(k1, 2) * besselI(k1, 0)/2 - besselI(k1, 1)^2
      der2 <- n * a/besselI(k1, 0)^2
      k2 <- k1 + der/der2
      while ( abs(k2 - k1) > 1e-08 ) {
        k1 <- k2
        der <- con - n * besselI(k1, 1, expon.scaled = TRUE)/besselI(k1, 0, expon.scaled = TRUE)
        a <- besselI(k1, 0)^2/2 + besselI(k1, 2) * besselI(k1, 0)/2 - besselI(k1, 1)^2
        der2 <- n * a/besselI(k1, 0)^2
        k2 <- k1 + der/der2
      }
    } else k2 <- k1
    kapa[i] <- k2
  }

  Rw <- sqrt( sum(kapa * Ri * cos(mi))^2 + sum( kapa * Ri * sin(mi))^2 )
  Ta <- 2 * (sum(kapa * Ri) - Rw)

  pta <- numeric(B)
  for (j in 1:B) {
    id <- sample(N, N)
    C <- Rfast::group(x1[id], ina)
    S <- Rfast::group(x2[id], ina)
    mi <- atan(S/C) + pi * as.numeric(C < 0)
    Ri <- sqrt(C^2 + S^2)
    ki <- (1.28 - 0.53 * Ri^2/ni^2) * tan(0.5 * pi * Ri/ni)
    coni <- ki
    x <- u[id]
    for (i in 1:2) {
      n <- ni[i]
      coni[i] <- sum( cos( x[ina == i] - mi[i] ) )
      con <- coni[i]
      k1 <- ki[i]
      if (k1 < 710) {
        der <- con - n * besselI(k1, 1, expon.scaled = TRUE)/besselI(k1, 0, expon.scaled = TRUE)
        a <- besselI(k1, 0)^2/2 + besselI(k1, 2) * besselI(k1, 0)/2 - besselI(k1, 1)^2
        der2 <- n * a/besselI(k1, 0)^2
        k2 <- k1 + der/der2
        while ( abs(k2 - k1) > 1e-08 ) {
          k1 <- k2
          der <- con - n * besselI(k1, 1, expon.scaled = TRUE)/besselI(k1, 0, expon.scaled = TRUE)
          a <- besselI(k1, 0)^2/2 + besselI(k1, 2) * besselI(k1, 0)/2 - besselI(k1, 1)^2
          der2 <- n * a/besselI(k1, 0)^2
          k2 <- k1 + der/der2
        }
      } else k2 <- k1
      kapa[i] <- k2
    }
    Rw <- sqrt( sum(kapa * Ri * cos(mi))^2 + sum( kapa * Ri * sin(mi))^2 )
    pta[j] <- 2 * (sum(kapa * Ri) - Rw)
  }

  pvalue <- ( sum(pta > Ta) + 1 ) / (B + 1) 
  res <- c(Ta, pvalue)
  names(res) <- c("test", "p-value")
  res
}