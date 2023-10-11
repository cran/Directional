hetcirc.perm <- function(u1, u2, rads = TRUE, B = 999) {

  if ( !rads )  {
    u1 <- u1 * pi/180
    u2 <- u2 * pi/180
  }
  u <- c(u1, u2)
  ina <- c( rep(1, length(u1) ), rep(2, length(u2) ) )
  N <- length(u)
  ni <- tabulate(ina)
  kapa <- numeric(2)
  x1 <- cos(u)   ;   x2 <- sin(u)
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
    id <- Rfast2::Sample.int(N, N)
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

  p.value <- ( sum(pta > Ta) + 1 ) / (B + 1)
  statistic <- Ta       ;   names(statistic) <- "het test statistic"
  parameter <- "NA"     ;   names(parameter) <- "df"
  alternative <- "The 2 circular means differ"
  method <- "Permutation ANOVA for 2 circular means using the heterogeneous approach"
  data.name <- c("data ", " groups")
  result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                  alternative = alternative, method = method, data.name = data.name )
  class(result) <- "htest"
  return(result)
}
