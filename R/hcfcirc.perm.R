hcfcirc.perm <- function(u1, u2, rads = TRUE, B = 999) {
  if ( !rads )  {
    u1 <- u1 * pi/180
    u2 <- u2 * pi/180
  }
  u <- c(u1, u2)
  ina <- c( rep(1, length(u1) ), rep(2, length(u2) ) )
  ni <- tabulate(ina)
  n <- sum(ni)
  x1 <- cos(u)
  x2 <- sin(u)
  Ci <- Rfast::group(x1, ina)
  Si <- Rfast::group(x2, ina)
  Ri <- sqrt(Ci^2 + Si^2)
  V <- sum(Ri)
  C <- sum(Ci)
  S <- sum(Si)
  R <- sqrt(C^2 + S^2)

  mu <- atan(S/C) + pi * ( C < 0 )
  con <- sum( cos(u - mu) )
  k1 <- (1.28 - 0.53 * R^2/n^2) * tan(0.5 * pi * R/n)
  if (k1 < 710) {
    der <- con - n * besselI(k1, 1, expon.scaled = TRUE)/besselI(k1, 0, expon.scaled = TRUE)
    a <- besselI(k1, 0)^2/2 + besselI(k1, 2) * besselI(k1, 0)/2 - besselI(k1, 1)^2
    der2 <- n * a/besselI(k1, 0)^2
    k2 <- k1 + der/der2
    while (abs(k1 - k2) > 1e-7) {
      k1 <- k2
      der <- con - n * besselI(k1, 1, expon.scaled = TRUE)/besselI(k1, 0, expon.scaled = TRUE)
      a <- besselI(k1, 0)^2/2 + besselI(k1, 2) * besselI(k1, 0)/2 - besselI(k1, 1)^2
      der2 <- n * a/besselI(k1, 0)^2
      k2 <- k1 + der/der2
    }
  } else k2 <- k1
  kapa <- k2

  if (kapa > 2) {
    Ft <- (n - 2) * (V - R)/(n - V)
  } else if (kapa < 2 & kapa > 1) {
    Ft <- (1 + 3/(8 * kapa) ) * (n - 2) * (V - R) / (n - V)
  } else  Ft <- NA

  pvalue <- NA
  if ( !is.na(Ft) ) {
    pft <- numeric(B)
    for (i in 1:B) {
      id <- sample(n, n)
      Ci <- Rfast::group(x1[id], ina)
      Si <- Rfast::group(x2[id], ina)
      Ri <- sqrt(Ci^2 + Si^2)
      V <- sum(Ri)

      if (kapa > 2) {
        pft[i] <- (n - 2) * (V - R)/(n - V)
      } else if (kapa < 2 & kapa > 1) {
        pft[i] <- (1 + 3/(8 * kapa)) * (n - 2) * (V - R) / (n - V)
      } else  pft[i] <- NA
    }
    pvalue <- ( sum(pft > Ft) + 1 ) / (B + 1)
  }
  res <- c(Ft, pvalue, kapa)
  names(res) <- c("test", "p-value", "kappa")
  res
}
