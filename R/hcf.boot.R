hcf.boot <- function(x1, x2, fc = TRUE, B = 999) {

  n1 <- dim(x1)[1]    ;  n2 <- dim(x2)[1]
  x <- rbind(x1, x2)
  ina <- c( rep(1, n1), rep(2, n2) )
  ni <- c(n1, n2)
  p <- dim(x)[2]  ## dimensionality of the data
  n <- n1 + n2  ## sample size of the data
  S <- rowsum(x, ina)
  m1 <- S[1, ]   ;   m2 <- S[2, ]
  m1 <- m1 / sqrt( sum(m1^2) )
  m2 <- m2 / sqrt( sum(m2^2) )
  Ri <- sqrt( Rfast::rowsums(S^2) )  ## the resultant length of each group
  S <- Rfast::colsums(x)
  R <- sqrt( sum(S^2) )  ## the resultant length based on all the data
  ## Next we stimate the common concentration parameter kappa
  kapaa <- Directional::vmf.mle(x, fast = TRUE)$kappa
  m <- S / R
  ## kapaa is the estimated concentration parameter based on all the data
  Ft <- (n - 2) * ( sum(Ri) - R) / ( n - sum(Ri) )
  if (fc) {  ## correction is used
    if (p == 3) {
      Ft <- kapaa * (1/kapaa - 1/(5 * kapaa^3)) * Ft
    } else if (p > 3)  Ft <- kapaa * ( 1/kapaa - (p - 3)/(4 * kapaa^2) - (p - 3)/(4 * kapaa^3) ) * Ft
  }

  rot1 <- t( Directional::rotation(m1, m) )
  rot2 <- t( Directional::rotation(m2, m) )
  y1 <- x1 %*% rot1
  y2 <- x2 %*% rot2
  ftb <- numeric(B)

  for (i in 1:B) {
    b1 <- sample(n1, n1, replace = TRUE)
    b2 <- sample(n2, n2, replace = TRUE)
    yb <- rbind(y1[b1, ], y2[b2, ])
    S <- rowsum(yb, ina)
    Ri <- sqrt( Rfast::rowsums(S^2) )
    S <- Rfast::colsums(yb)
    R <- sqrt( sum(S^2) )
    kapa <- Directional::vmf.mle(x, fast = TRUE)$kappa
    ftb[i] <- (n - 2) * ( sum(Ri) - R) / ( n - sum(Ri) )
    if (fc) {  ## correction is used
      if (p == 3) {
        ftb[i] <- kapa * (1/kapa - 1/(5 * kapa^3)) * ftb[i]
      } else if (p > 3)  ftb[i] <- kapa * ( 1/kapa - (p - 3)/(4 * kappa^2) - (p - 3)/(4 * kapa^3) ) * ftb[i]
    }
  }

  pvalue <- ( sum(ftb > Ft) + 1 ) / (B + 1)
  res <- c(Ft, pvalue, kapaa)
  names(res) <- c('test', 'p-value', 'kappa')
  res
}
