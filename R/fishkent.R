################################
#### H0: Fisher versus H1: Kent distribution
#### Tsagris Michail 06/2014
#### mtsagris@yahoo.gr
#### References: Louis-Paul Rivest (1986)
#### Statistics & Probability Letters, 4: 1-4.
#### Modified Kent's statistics for testing goodness of fit for the
#### Fisher distribution in small concentrated samples
################################
fishkent <- function(x, B = 999) {
  ## x contains the data
  ## B is by default eaual to 999 bootstrap re-samples
  ## If B==1 then no bootstrap is performed
  n <- dim(x)[1]  ## sample size
  estim <- Directional::vmf.mle(x, fast = TRUE)
  k <- estim$kappa  ## the estimated concentration parameter
  ## under the H0, that the Fisher distribution is true
  mu <- estim$mu  ## the estimated mean direction under H0
  e1 <- c(1, 0, 0)
  i3 <- diag(3)
  P <- i3 -  tcrossprod(e1 - mu) / (1 - mu[1])
  y <- tcrossprod(x, P)[, 2:3]
  lam <- eigen( crossprod(y) )$values / n
  rat <- besselI(k, 0.5, expon.scaled = TRUE) / besselI(k, 2.5, expon.scaled = TRUE)
  Ta <- n * (k / 2)^2 * rat * (lam[1] - lam[2])^2

  if (B == 1) {
    p.value <- pchisq(Ta, 2, lower.tail = FALSE)
    parameter <- 2     ;   names(parameter) <- "df"
    statistic <- Ta  ;   names(statistic) <- "Test statistic"
    alternative <- "Kent is prefered to von Mises-Fisher"
    method <- "Asymptotic rotational symmetry (von Mises-Fisher versus Kent distribution)"
    data.name <- c("data")
    result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                    alternative = alternative, method = method, data.name = data.name )
    class(result) <- "htest"

  } else {
    Tb <- numeric(B)
    for (i in 1:B) {
      nu <- Rfast2::Sample.int(n, n, replace = TRUE)
      z <- x[nu, ]
      estim <- Directional::vmf.mle(z, fast = TRUE)
      k <- estim$kappa  ## the estimated concentration parameter
      ## under the H0, that the Fisher distribution is true
      mu <- estim$mu  ## the estimated mean direction under H0
      P <- i3 -  tcrossprod(e1 - mu) / (1 - mu[1])
      y <- tcrossprod(z, P)[, 2:3]
      lam <- eigen( crossprod(y) )$values/n
      rat <- besselI(k, 0.5, expon.scaled = TRUE) / besselI(k, 2.5, expon.scaled = TRUE)
      Tb[i] <- n * ( k / 2 )^2 * rat * (lam[1] - lam[2])^2
    }
    p.value <-  (sum(Tb > Ta) + 1) / (B + 1)
    parameter <- "NA"     ;   names(parameter) <- "df"
    statistic <- Ta  ;   names(statistic) <- "Test statistic"
    alternative <- "Kent is prefered to von Mises-Fisher"
    method <- "Bootstrap rotational symmetry (von Mises-Fisher versus Kent distribution)"
    data.name <- c("data")
    result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                    alternative = alternative, method = method, data.name = data.name )
    class(result) <- "htest"
  }

  return(result)
}
