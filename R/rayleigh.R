################################
#### Rayleigh test of uniformity
#### Tsagris Michail 6/2014
#### mtsagris@yahoo.gr
#### References: Mardia K.V., Kent J.T. and Bibby J.M. (1979) pg 439-440.  Multivariate analaysis
#### Mardia Kanti V. and Jupp Peter E. (2000) pg. 94-95. Directional statistics
################################

rayleigh <- function(x, modif = TRUE, B = 999) {
  ## x contains the data in Euclidean coordinates
  ## B is by default eaual to 999 bootstrap samples
  ## If B==1 then no bootstrap is performed

  x <- as.matrix(x)  ## makes sure x is a matrix
  x <- x/sqrt(rowSums(x^2))  ## makes sure x contains unit vectors
  p <- ncol(x)  ## dimensionality of the data
  n <- nrow(x)  ## sample size of the data
  m <- colSums(x)
  T <- sum( m^2 ) *p / n

  if (modif == TRUE) {
    T <- ( 1 - 1/(2 * n) ) * T + T^2 / ( 2 * n * (p + 2) )
  }

  if (B == 1) {
    pvalue <- 1 - pchisq(T, p)
    res <- c(T, pvalue)
    names(res) <- c('test', 'p-value')

  } else {
    Tb <- numeric(B)
    for (i in 1:B) {
      x <- matrix( rnorm(p * n), ncol = p )
      x <- x / sqrt( rowSums(x^2) )
      mb <- colSums(x)
      Tb[i] <- p * sum( mb^2 ) / n
    }

    res <- c( T, (sum(Tb > T) + 1)/(B + 1) )
    names(res) <- c('test', 'Bootstrap p-value')
  }

  res

}
