het.perm <- function(x1, x2, B = 999) {

  n1 <- dim(x1)[1]    ;  n2 <- dim(x2)[1]
  x <- rbind(x1, x2)
  ina <- c( rep(1, n1), rep(2, n2) )
  ni <- c(n1, n2)
  p <- dim(x)[2]  ## dimensionality of the data
  n <- n1 + n2  ## sample size of the data
  S <- rowsum(x, ina)
  kapa <- numeric(2)
  mi <- rowsum(x, ina) / ni
  kapa[1] <- Directional::vmf( x1, fast = TRUE )$kappa
  kapa[2] <- Directional::vmf( x2, fast = TRUE )$kappa
  tw <- Rfast::colsums(kapa * ni * mi)
  Tt <- sum( kapa * ni * sqrt( Rfast::rowsums(mi^2) ) ) - sqrt( sum(tw^2) )

  ptt <- numeric(B)
  for (i in 1:B) {
    ind <- sample(ina, n)
    S <- rowsum(x, ind)
    kapa <- numeric(2)
    mi <- rowsum(x, ind) / ni
    kapa[1] <- Directional::vmf( x[ind == 1, ], fast = TRUE )$kappa
    kapa[2] <- Directional::vmf( x[ind == 2, ], fast = TRUE )$kappa
    tw <- Rfast::colsums(kapa * ni * mi)
    ptt[i] <- sum( kapa * ni * sqrt( Rfast::rowsums(mi^2) ) ) - sqrt( sum(tw^2) )
  }
  pvalue <- ( sum(ptt > Tt) + 1 ) / (B + 1)
  res <- c(2 * Tt, pvalue)
  names(res) <- c('test', 'p-value')
  res
}
