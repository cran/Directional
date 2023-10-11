################################
#### Hypothesis test for a mean direction
#### Tsagris Michail 6/2014
#### mtsagris@yahoo.gr
#### References: Mardia Kanti V. and Jupp Peter E. (2000)
#### Directional statistics, pg. 212
################################
meandir.test <- function(x, mu, B = 999) {
  ## x is the sample
  ## mu is the hypothesized mean direction under H0
  p <- dim(x)[2]  ## dimensionality of the data
  n <- dim(x)[1]  ## sample size of the data
  k1 <- Directional::vmf.mle(x)$k  ## concentration parameter under H1
  xbar <- Rfast::colmeans(x)  ## x-bar
  m1 <- xbar / sqrt( sum(xbar^2) )
  sxm <- sum(x %*% mu)

  lik <- function(k, x, sxm)  n * (p/2 - 1) * log(k) + k * sxm -  n * ( log( besselI(k, p/2 - 1, expon.scaled = TRUE) ) + k )

  qa0 <- optimize(lik, c(0, 100000), x = x, sxm = sxm, maximum = TRUE)
  k0 <- qa0$maximum  ## concentration parameter under H0
  apk0 <- (1 - p/2) * log(k0/2) + lgamma(p/2) + log( besselI(k0, p/2 - 1, expon.scaled = TRUE) ) + k0
  apk1 <- (1 - p/2) * log(k1/2) + lgamma(p/2) + log( besselI(k1, p/2 - 1, expon.scaled = TRUE) ) + k1

  w <- 2 * n * ( k1 * sqrt( sum(xbar^2) ) - k0 * sum(mu * xbar) - apk1 + apk0 )
  p.value <- pchisq(w, p - 1, lower.tail = FALSE)
  parameter <- p - 1     ;   names(parameter) <- "df"
  statistic <- w  ;   names(statistic) <- "Test statistic"
  alternative <- "Mean direction is not equal to some predefined direction"
  method <- "Hypothesis test for a mean direction"
  data.name <- c("data")
  result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                  alternative = alternative, method = method, data.name = data.name )
  class(result) <- "htest"

  if (B > 1) {
    A <- Directional::rotation(m1, mu)
    y <- tcrossprod(x, A)  ## bring the data under H0
    ## y has mean direction equal to mu
    wb <- numeric(B)

    for (i in 1:B) {
      nu <- Rfast2::Sample.int(n, n, replace = TRUE)
      z <- y[nu, ]
      k1 <- Directional::vmf.mle(z, fast = TRUE)$k  ## concentration parameter under H1
      zbar <- Rfast::colmeans(z)  ## z-bar
      sxm <- sum(z %*% mu)
      qa0 <- optimize(lik, c(0, 100000), x = z, sxm = sxm, maximum = TRUE)
      k0 <- qa0$maximum  ## concentration parameter under H0
      apk0 <- (1 - p/2) * log(k0/2) + lgamma(p/2) + log( besselI(k0, p/2 - 1, expon.scaled = TRUE) ) + k0
      apk1 <- (1 - p/2) * log(k1/2) + lgamma(p/2) + log( besselI(k1, p/2 - 1, expon.scaled = TRUE) ) + k1
      wb[i] <- 2 * n * ( k1 * sqrt( sum(zbar^2) ) - k0 * sum(mu * zbar) - apk1 + apk0 )
    }

    p.value <- (sum(wb > w) + 1)/(B + 1)
    parameter <- "NA"     ;   names(parameter) <- "df"
    statistic <- w  ;   names(statistic) <- "Test statistic"
    alternative <- "Mean direction is not equal to some predefined direction"
    method <- "Bootstrap hypothesis test for a mean direction"
    data.name <- c("data")
    result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                    alternative = alternative, method = method, data.name = data.name )
    class(result) <- "htest"
  }

  return(result)
}
