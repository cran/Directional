################################
#### Watson test of uniformity with circular data
#### Tsagris Michail 3/2016
#### mtsagris@yahoo.gr
#### References: Jammalamadaka, S. Rao and SenGupta, A. (2001).
#### Topics in Circular Statistics, pg. 156-157
################################
watson <- function(u, rads = FALSE, R = 1) {
  ## u is a vector with circular data
  ## if data are in rads set it to TRUE
  ## R is for Monte Carlo estimate of the p-value
  if ( !rads )   u <- u / 180 * pi
  u <- Rfast::Sort(u) / (2 * pi)
  n <- length(u)
  i <- 1:n
  Wn <- sum( ( ( u - (i - 0.5)/n ) - ( sum(u) / n - 0.5 ) )^2 ) + 1 / ( 12 * n )

  if ( R == 1 ) {  ## asymptotic p-value is returned
    m <- 1:20
    p.value <- 2 * sum( ( - 1 )^( m - 1 ) * exp(-2 * m^2 * pi^2 * Wn) )
  } else {
    x <- matrix( Rfast2::Runif(n * R, 0, 2 * pi), ncol = R )
    x <- Rfast::colSort(x) / (2 * pi)
    bwn <- Rfast::colMaxs(x - (i - 1)/n, value = TRUE) + Rfast::colMaxs(i/n - x, value = TRUE)
    p.value <- ( sum(bwn > Wn) + 1 ) / (R + 1)
  }

  parameter <- "NA"     ;   names(parameter) <- "df"
  statistic <- Wn  ;   names(statistic) <- "Test statistic"
  alternative <- "The distribution is circular uniform"
  method <- "Watson test of uniformity with circular data"
  data.name <- c("data")
  result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                  alternative = alternative, method = method, data.name = data.name )
  class(result) <- "htest"
}
