embedcirc.perm <- function(u1, u2, rads = TRUE, B = 999) {

  u <- c(u1, u2)
  if ( !rads )  u <- u * pi/180
  ina <- c( rep(1, length(u1) ), rep(2, length(u2) ) )
  ni <- tabulate(ina)
  n <- sum(ni)
  g <- 2
  x1 <- cos(u)
  x2 <- sin(u)
  Ci <- Rfast::group(x1, ina)
  Si <- Rfast::group(x2, ina)
  Rbi <- (Ci^2 + Si^2)/ni^2
  C <- sum(Ci)
  S <- sum(Si)
  Rbar <- sqrt(C^2 + S^2)/n

  mu <- atan(S/C) + pi * (C < 0)
  con <- sum( cos(u - mu) )
  k1 <- (1.28 - 0.53 * Rbar^2) * tan(0.5 * pi * Rbar)
  if (k1 < 710) {
    der <- con - n * besselI(k1, 1, expon.scaled = TRUE)/besselI(k1, 0, expon.scaled = TRUE)
    a <- besselI(k1, 0)^2/2 + besselI(k1, 2) * besselI(k1, 0)/2 - besselI(k1, 1)^2
    der2 <- n * a/besselI(k1, 0)^2
    k2 <- k1 + der/der2
    while (abs(k1 - k2) > 1e-07) {
      k1 <- k2
      der <- con - n * besselI(k1, 1, expon.scaled = TRUE)/besselI(k1, 0, expon.scaled = TRUE)
      a <- besselI(k1, 0)^2/2 + besselI(k1, 2) * besselI(k1, 0)/2 - besselI(k1, 1)^2
      der2 <- n * a/besselI(k1, 0)^2
      k2 <- k1 + der/der2
    }
  } else k2 <- k1
  kapa <- k2

  Fb <- (sum(ni * Rbi) - n * Rbar^2) / ( (n - sum(ni * Rbi) ) / (n - 2) )
  Fc <- (1 - 1/(5 * kapa) - 1/( 10 * kapa^2) ) * Fb

  pfb <- numeric(B)
  for (i in 1:B) {
    ind <- sample(ina, n)
    Ci <- Rfast::group(x1, ind)
    Si <- Rfast::group(x2, ind)
    Rbi <- (Ci^2 + Si^2)/ni^2
    pfb[i] <- (sum(ni * Rbi) - n * Rbar^2) / ( (n - sum(ni * Rbi) ) / (n - 2) )
  }
  pvalue <- ( sum(pfb > Fb) + 1 ) / (B + 1)
  res <- c(Fc, pvalue, kapa)
  names(res) <- c("test", "p-value", "kappa")
  res
}
