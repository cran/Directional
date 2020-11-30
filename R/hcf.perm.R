hcf.perm <- function(x1, x2, B = 999) {

  n1 <- dim(x1)[1]    ;  n2 <- dim(x2)[1]
  x <- rbind(x1, x2)
  ina <- c( rep(1, n1), rep(2, n2) )
  ni <- c(n1, n2)
  p <- dim(x)[2]  ## dimensionality of the data
  n <- n1 + n2  ## sample size of the data
  S <- rowsum(x, ina)
  Ri <- sqrt( Rfast::rowsums(S^2) )  ## the resultant length of each group
  S <- Rfast::colsums(x)
  R <- sqrt( sum(S^2) )  ## the resultant length based on all the data

  Ft <- (n - 2) * (sum(Ri) - R) / ( n - sum(Ri) )

  pft <- numeric(B)
  for (i in 1:B) {
    ind <- sample(ina, n)
    S <- rowsum(x, ind)
    Ri <- sqrt( Rfast::rowsums(S^2) )
    pft[i] <- (n - 2) * (sum(Ri) - R) / ( n - sum(Ri) )
  }

  pvalue <- ( sum(pft > Ft) + 1 ) / (B + 1)
  res <- c(Ft, pvalue)
  names(res) <- c('test', 'p-value')
  res
}
