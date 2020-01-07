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
	  if ( identical(class(x)[1], "FBM") ) {
	    colMeans_sub <- function(X, ind) colMeans(X[, ind])
	    m1 <- bigstatsr::big_apply(x, a.FUN = colMeans_sub, a.combine = 'c')
    } else  m1 <- Rfast::colsums(x)
    R <- sqrt(sum(m1^2))/n
    m <- m1/n/R
    k <- R * (p - R^2)/(1 - R^2)
    if (k < 1e+05) {
        lik1 <- (p/2 - 1) * log(k) - log(besselI(k, p/2 - 1,
            expon.scaled = TRUE)) - k + k * R
        apk <- Apk(p, k)
        k <- k - (apk - R)/(1 - apk^2 - (p - 1)/k * apk)
        lik2 <- (p/2 - 1) * log(k) - log(besselI(k, p/2 - 1,
            expon.scaled = TRUE)) - k + k * R
        while (lik2 - lik1 > tol) {
            lik1 <- lik2
            apk <- Apk(p, k)
            k <- k - (apk - R)/(1 - apk^2 - (p - 1)/k * apk)
            lik2 <- (p/2 - 1) * log(k) - log(besselI(k, p/2 -
                1, expon.scaled = TRUE)) - k + k * R
        }
    }
    else k <- k
    loglik <- n * lik2 - 0.5 * n * p * log(2 * pi)
    vark <- 1 / ( n * (1 - Apk(p, k)/k - Apk(p, k)^2) )
    res <- list(mu = m, kappa = k, MRL = R, vark = vark, loglik = loglik)
  }

  res
}
