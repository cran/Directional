embedcirc.perm <- function(u1, u2, rads = TRUE, B = 999) {

  if ( !rads )  {
    u1 <- u1 * pi/180
    u2 <- u2 * pi/180
  }
  u <- c(u1, u2)

  ina <- c( rep(1, length(u1) ), rep(2, length(u2) ) )
  ni <- tabulate(ina)
  n <- sum(ni)
  g <- 2
  x1 <- cos(u)   ;    x2 <- sin(u)
  Ci <- Rfast::group(x1, ina)
  Si <- Rfast::group(x2, ina)
  Rbi <- (Ci^2 + Si^2)/ni^2
  C <- sum(Ci)
  S <- sum(Si)
  Rbar <- sqrt(C^2 + S^2)/n

  Fb <- (n - 2) * ( sum(ni * Rbi) - n * Rbar^2 ) / ( n - sum(ni * Rbi) )

  pfb <- numeric(B)
  for (i in 1:B) {
    ind <- sample(ina, n)
    Ci <- Rfast::group(x1, ind)
    Si <- Rfast::group(x2, ind)
    Rbi <- (Ci^2 + Si^2)/ni^2
    C <- sum(Ci)
    S <- sum(Si)
    Rbar <- sqrt(C^2 + S^2)/n
    pfb[i] <- (n - 2) * (sum(ni * Rbi) - n * Rbar^2) / ( n - sum(ni * Rbi) )
  }
  pvalue <- ( sum(pfb > Fb) + 1 ) / (B + 1)
  res <- c(Fb, pvalue)
  names(res) <- c("test", "p-value")
  res
}
