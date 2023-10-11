################################
#### Kuiper test of uniformity with circular data
#### Tsagris Michail 3/2016
#### mtsagris@yahoo.gr
#### References: Jammalamadaka, S. Rao and SenGupta, A. (2001).
#### Topics in Circular Statistics, pg. 153-155
################################
kuiper <- function(u, rads = FALSE, R = 1) {
  ## u is a vector with circular data
  ## if data are in rads set it to TRUE
  ## R is for Monte Carlo estimate of the p-value
  if ( !rads )   u <- u / 180 * pi
  u <- Rfast::Sort(u) / (2 * pi)
  n <- length(u)
  i <- 1:n
  f <- sqrt(n)
  Vn <- f * ( max(u - (i - 1)/n ) + max(i/n - u) )

  if ( R == 1 ) {  ## asymptotic p-value is returned
    m <- (1:50)^2
    a1 <- 4 * m * Vn^2
    a2 <- exp( -2 * m * Vn^2 )
    b1 <- 2 * ( a1 - 1 ) * a2
    b2 <- 8 * Vn / ( 3 * f ) * m * (a1 - 3) * a2
    p.value <- sum(b1 - b2)
  } else {
    x <- matrix( Rfast2::Runif(n * R, 0, 2 * pi), ncol = R)
    x <- Rfast::colSort(x) / (2 * pi)
    bvn <- f * ( Rfast::colMaxs(x - (i - 1)/n, value = TRUE) + Rfast::colMaxs(i/n - x, value = TRUE) )
    p.value <- ( sum(bvn > Vn) + 1 ) / (R + 1)
  }

  parameter <- "NA"     ;   names(parameter) <- "df"
  statistic <- Vn  ;   names(statistic) <- "Test statistic"
  alternative <- "The distribution is circular uniform"
  method <- "Kuiper test of uniformity with circular data"
  data.name <- c("data")
  result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                  alternative = alternative, method = method, data.name = data.name )
  class(result) <- "htest"
}





