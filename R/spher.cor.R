################################
#### Spherical-spherical correlation
#### Tsagris Michail 4/2014
#### mtsagris@yahoo.gr
#### References: Kanti V. Mardia and Peter E. Jupp
#### Directional statistics p.g. 254-255
################################

spher.cor <- function(x, y) {
  ## x and y are two (hyper-)spherical variables
  
  x <- as.matrix(x)
  y <- as.matrix(y)
  x <- x / sqrt( Rfast::rowsums(x^2) )
  y <- y / sqrt( Rfast::rowsums(y^2) )
  
  p <- dim(x)[2]  ## dimension of x
  q <- dim(y)[2]  ## dimension of y
  n <- dim(x)[1]  ## sample size
  x <- t(x) - Rfast::colmeans(x)     ## subtract the mean
  y <- t(y) - Rfast::colmeans(y)   ## subtract the mean

  s11 <- tcrossprod(x) / n
  s12 <- tcrossprod( x, y ) / n
  s21 <- t( s12 )
  s22 <- tcrossprod(y) / n
  
  a1 <- solve(s11, s12)
  a2 <- solve(s22, s21) 
  rsq <- sum( t(a1) * a2)
  test <- n * rsq
  pvalue <- pchisq(test, p * q, lower.tail = FALSE)
  res <- c(rsq, pvalue)
  
  names(res) <- c('R-squared', 'p-value')
  res
  
}