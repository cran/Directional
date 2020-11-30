embed.perm <- function(x1, x2, B = 999) {
  n1 <- dim(x1)[1]    ;  n2 <- dim(x2)[1]
  x <- rbind(x1, x2)
  ina <- c( rep(1, n1), rep(2, n2) )
  ni <- c(n1, n2)
  p <- dim(x)[2]  ## dimensionality of the data
  n <- n1 + n2  ## sample size of the data
  S <- rowsum(x, ina) / ni
  Rbi <- sqrt( Rfast::rowsums(S^2) ) ## the mean resultant length of each group
  S <- Rfast::colmeans(x)
  Rbar <- sqrt( sum(S^2) )  ## the mean resultant length based on all the data
  Ft <- (n - 2) * (p - 1) * ( sum(ni * Rbi^2) - n * Rbar^2) / ( (p - 1) * ( n - sum(ni * Rbi^2) ) )

  pft <- numeric(B)
  for (i in 1:B) {
    ind <- sample(ina, n)
    S <- rowsum(x, ind) / ni
    Rbi <- sqrt( Rfast::rowsums(S^2) )
    pft[i] <- (n - 2) * ( sum(ni * Rbi^2) - n * Rbar^2) / ( n - sum(ni * Rbi^2) )
  }

  pvalue <- ( sum(pft > Ft) + 1 ) / (B + 1)
  res <- c(Ft, pvalue)
  names(res) <- c("test", "p-value")
  res
}
