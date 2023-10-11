lr.aov <- function(x, ina) {

  ina <- as.numeric(ina)
  g <- max(ina)  ## how many groups are there
  p <- dim(x)[2]  ## dimensionality of the data
  n <- dim(x)[1]  ## sample size of the data
  S <- rowsum(x, ina)
  Ri <- sqrt( Rfast::rowsums(S^2) )  ## the resultant length of each group
  S <- Rfast::colsums(S)
  R <- sqrt( sum(S^2) )  ## the resultant length based on all the data

  Apk <- function(p, k)  besselI(k, p/2, expon.scaled = TRUE) / besselI(k, p/2 - 1, expon.scaled = TRUE)

  ## Next we stimate the common concentration parameter kappa under H0 and H1
  Rk <- R/n
  k1 <- Rk * (p - Rk^2)/(1 - Rk^2)
  k2 <- k1 - (Apk(p, k1) - Rk) / ( 1 - Apk(p, k1)^2 - (p - 1)/k1 * Apk(p, k1) )
  while ( abs(k2 - k1) > 1e-07 ) {
    k1 <- k2
    k2 <- k1 - (Apk(p, k1) - Rk) / (1 - Apk(p, k1)^2 - (p - 1)/k1 * Apk(p, k1) )
  }
  k0 <- k2  ## concentration parameter under H0

  Rk <- sum(Ri)/n
  k1 <- Rk * (p - Rk^2)/(1 - Rk^2)
  k2 <- k1 - (Apk(p, k1) - Rk) / ( 1 - Apk(p, k1)^2 - (p - 1)/k1 * Apk(p, k1) )

  while ( abs(k2 - k1) > 1e-07 ) {
    k1 <- k2
    k2 <- k1 - (Apk(p, k1) - Rk) / ( 1 - Apk(p, k1)^2 - (p - 1)/k1 * Apk(p, k1) )
  }
  k1 <- k2  ## concentration parameter under H1

  apk0 <- (1 - p/2) * log(k0/2) + lgamma(p/2) + log( besselI( k0, p/2 - 1, expon.scaled = TRUE ) ) + k0
  apk1 <- (1 - p/2) * log(k1/2) + lgamma(p/2) + log( besselI( k1, p/2 - 1, expon.scaled = TRUE ) ) + k1
  w <- 2 * (k1 * sum(Ri) - k0 * R - n * apk1 + n * apk0)
  p.value <- pchisq(w, (g - 1) * (p - 1), lower.tail = FALSE)

  statistic <- w   ;   names(statistic) <- "chi-square test statistic"
  parameter <- (g - 1) * (p - 1)     ;   names(parameter) <- "df"
  alternative <- "At least one directional mean vector differs"
  method <- "ANOVA for directional data using the log-likelihood ratio test"
  data.name <- c("data ", " groups")
  result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                  alternative = alternative, method = method, data.name = data.name )
  class(result) <- "htest"
  return(result)
}

