################################
#### Spherical-spherical correlation
#### Tsagris Michail 4/2014
#### mtsagris@yahoo.gr
#### References: Kanti V. Mardia and Peter E. Jupp
#### Directional statistics p.g. 254-255
################################
spher.cor <- function(x, y) {
  ## x and y are two (hyper-)spherical variables
  p <- dim(x)[2]  ## dimension of x
  q <- dim(y)[2]  ## dimension of y
  n <- dim(x)[1]  ## sample size
  x <- Rfast::eachrow( x, Rfast::colmeans(x), oper = "-" )     ## subtract the mean
  y <- Rfast::eachrow(y, Rfast::colmeans(y), oper = "-" )   ## subtract the mean
  s11 <- crossprod(x) / n
  s12 <- crossprod( x, y ) / n
  s21 <- t( s12 )
  s22 <- crossprod(y) / n
  a1 <- solve(s11, s12)
  a2 <- solve(s22, s21)
  rsq <- sum( t(a1) * a2)
  test <- n * rsq
  pvalue <- pchisq(test, p * q, lower.tail = FALSE)
  res <- c(rsq, pvalue)
  names(res) <- c('R-squared', 'p-value')
  res
}
