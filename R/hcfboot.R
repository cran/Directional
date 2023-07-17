hcfboot <- function(x, ina, B = 999) {
  
  x <- x[order(ina), ]
  ina <- as.numeric(ina)
  ina <- sort(ina)
  ni <- tabulate(ina)
  k <- length(ni)
  dm <- dim(x)
  p <- dm[2]  ## dimensionality of the data
  n <- dm[1]  ## sample size of the data
  S <- rowsum(x, ina)
  Ri <- sqrt( Rfast::rowsums(S^2) )  ## the resultant length of each group
  mi <- S / Ri
  S <- Rfast::colsums(x)
  R <- sqrt( sum(S^2) )  ## the resultant length based on all the data
  ## Next we stimate the common concentration parameter kappa
  m <- S / R
  kapaa <- Directional::vmf.mle(x, fast = TRUE)$kappa
  ## kapaa is the estimated concentration parameter based on all the data
  Ft <- (n - k) / (k - 1) * ( sum(Ri) - R) / ( n - sum(Ri) )
  
  y <- list()
  for (j in 1:k) {
    rot <- t( Directional::rotation(mi[j, ], m) )
    y[[ j ]] <- x[ina == j, ] %*% rot 
  }
  
  ftb <- numeric(B)

  for (i in 1:B) {
    yb <- NULL
    for (j in 1:k) {
      b <- Rfast2::Sample.int(ni[j], ni[j], replace = TRUE)
      yb <- rbind( yb, y[[ j ]][b, ] )
    }
    S <- rowsum(yb, ina)
    Ri <- sqrt( Rfast::rowsums(S^2) )
    S <- Rfast::colsums(yb)
    R <- sqrt( sum(S^2) )
    ftb[i] <- (n - k) / (k - 1) * ( sum(Ri) - R) / ( n - sum(Ri) )
  }

  pvalue <- ( sum(ftb > Ft) + 1 ) / (B + 1)
  res <- c(Ft, pvalue, kapaa)
  names(res) <- c('test', 'p-value', 'kappa')
  res
}






