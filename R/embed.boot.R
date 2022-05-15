embed.boot <- function(x1, x2, B = 999) {
  ## x contains all the data
  ## ina is an indicator variable of each sample
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
  Ft <- (n - 2) * ( sum(ni * Rbi^2) - n * Rbar^2) / ( n - sum(ni * Rbi^2) )

  m1 <- S[1, ]   ;   m2 <- S[2, ]
  m1 <- m1 / sqrt( sum(m1^2) )
  m2 <- m2 / sqrt( sum(m2^2) )
  m <- S / Rbar
  rot1 <- t( Directional::rotation(m1, m) )
  rot2 <- t( Directional::rotation(m2, m) )
  y1 <- x1 %*% rot1
  y2 <- x2 %*% rot2
  ftb <- numeric(B)

  for (i in 1:B) {
    b1 <- Rfast2::Sample.int(n1, n1, replace = TRUE)
    b2 <- Rfast2::Sample.int(n2, n2, replace = TRUE)
    yb <- rbind(y1[b1, ], y2[b2, ])
    S <- rowsum(yb, ina) / ni
    Rbi <- sqrt( Rfast::rowsums(S^2) ) ## the mean resultant length of each group
    S <- Rfast::colmeans(yb)
    Rbar <- sqrt( sum(S^2) )  ## the mean resultant length based on all the data
    ftb[i] <-  (n - 2) * ( sum(ni * Rbi^2) - n * Rbar^2) / ( n - sum(ni * Rbi^2) )
  }
  pvalue <- ( sum(ftb > Ft) + 1 ) / (B + 1)
  res <- c(Ft, pvalue)
  names(res) <- c('test', 'p-value')
  res
}
