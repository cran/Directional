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

  if (B > 1) {
    A <- rotation(m1, mu)
    y <- tcrossprod(x, A)  ## bring the data under H0
    ## y has mean direction equal to mu
    wb <- numeric(B)
    for (i in 1:B) {
      nu <- sample(1:n, n, replace = TRUE)
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
    pvalue <- (sum(wb > w) + 1)/(B + 1)
  } else  pvalue <- pchisq(w, p - 1, lower.tail = FALSE)

  list(mean.dir = m1, pvalue = pvalue)
}
