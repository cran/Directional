################################
#### MLE for a von-Mises Fisher distribution
#### Tsagris Michail 10/2013
#### mtsagris@yahoo.gr
#### References: Suvrit Sra (2012, Computational Statistics)
#### A short note on parameter approximation for von Mises-Fisher distributions:
#### and a fast implementation of Is(x)
################################

vmf <- function(x, fast = FALSE, tol = 1e-07) {
  ## x contains the data
  ## tol specifies the tolerance value for convergence
  ## when estimating the concentration parameter
  if ( fast ) {
    res <- Rfast::vmf.mle(x, tol = tol)
  } else {
    p <- dim(x)[2]  ## dimensionality of the data
    n <- dim(x)[1]  ## sample size of the data
    Apk <- function(p, k)  besselI(k, p/2, expon.scaled = TRUE) / besselI(k, p/2 - 1, expon.scaled = TRUE)
    m1 <- Rfast::colsums(x)
    R <- sqrt( sum(m1^2) )/n  ## mean resultant length
    m <- m1 / n / R
    k <- numeric(10)
    i <- 1
    k[i] <- R * (p - R^2)/(1 - R^2)
    if (k[i] > 100000) {
      k <- k[i]
    } else {
      i <- 2
      apk <- Apk(p, k[i - 1])
      k[i] <- k[i - 1] - (apk - R)/( 1 - apk^2 - (p - 1)/k[i - 1] * apk )
      while ( abs(k[i] - k[i - 1]) > tol ) {
        i <- i + 1
        apk <- Apk(p, k[i - 1])
        k[i] <- k[i - 1] - (apk - R)/( 1 - apk^2 - (p - 1)/k[i - 1] * apk )
      }
      k <- k[i]
    }
    loglik <- n * (p/2 - 1) * log(k) - 0.5 * n * p * log(2 * pi) - n * ( log( besselI(k, p/2 - 1, expon.scaled = TRUE) ) + k ) + k * sum(x %*% m)
    vark <- 1 / ( n * (1 - Apk(p, k)/k - Apk(p, k)^2) )
    res <- list(mu = m, kappa = k, MRL = R, vark = vark, loglik = loglik)
  }
  res
}
