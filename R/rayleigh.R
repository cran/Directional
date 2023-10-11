################################
#### Rayleigh test of uniformity
#### Tsagris Michail 6/2014
#### mtsagris@yahoo.gr
#### References: Mardia K.V., Kent J.T. and Bibby J.M. (1979) pg 439-440.  Multivariate analaysis
#### Mardia Kanti V. and Jupp Peter E. (2000) pg. 94-95. Directional statistics
################################
rayleigh <- function(x, modif = TRUE, B = 999) {
  ## x contains the data in Euclidean coordinates
  ## B is by default equal to 999 bootstrap samples
  ## If B==1 then no bootstrap is performed
  p <- dim(x)[2]  ## dimensionality of the data
  n <- dim(x)[1]  ## sample size of the data
  m <- Rfast::colsums(x)
  test <- sum( m^2 ) * p / n
  if ( modif )  test <- ( 1 - 1/(2 * n) ) * test + test^2 / ( 2 * n * (p + 2) )

  if (B == 1) {
    p.value <- pchisq(test, p, lower.tail = FALSE)
    parameter <- p     ;   names(parameter) <- "df"
    statistic <- test  ;   names(statistic) <- "Test statistic"
    alternative <- "The distribution is not uniform"
    method <- "Rayleigh test of uniformity"
    data.name <- c("data")
    result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                    alternative = alternative, method = method, data.name = data.name )
    class(result) <- "htest"

  } else {
    tb <- numeric(B)
    for (i in 1:B) {
      x <- Rfast::matrnorm(n, p)
      x <- x / sqrt( Rfast::rowsums(x^2) )
      mb <- Rfast::colsums(x)
      tb[i] <- p * sum( mb^2 ) / n
    }
    p.value <- ( sum(tb > test) + 1 ) / (B + 1)
    parameter <- "NA"     ;   names(parameter) <- "df"
    statistic <- test  ;   names(statistic) <- "Test statistic"
    alternative <- "The distribution is not uniform"
    method <- "Bootstrap Rayleigh test of uniformity"
    data.name <- c("data")
    result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                    alternative = alternative, method = method, data.name = data.name )
    class(result) <- "htest"
  }

  result
}
