ptest <- function(x, B = 100) {
  d <- dim(x)[2]
  u <- Directional::riag(B, numeric(d))
  z <- tcrossprod(x, u) 
  p <- numeric(B)
  for (i in 1:B)  
  p[i] <- ks.test(z[, i], "punif", -1, 1)$p.value
  p <- sort(p)
  i <- 1:B
  pval <- min(B / i * p )
  list(pvalues = p, pvalue = pval)
}
