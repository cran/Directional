hclr.perm <- function(x1, x2, B = 999) {

  n1 <- dim(x1)[1]    ;  n2 <- dim(x2)[1]
  x <- rbind(x1, x2)
  ina <- c( rep(1, n1), rep(2, n2) )
  ni <- c(n1, n2)
  p <- dim(x)[2]
  n <- n1 + n2
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

  I0 <- besselI( k0, p/2 - 1 )
  I1 <- besselI( k1, p/2 - 1 )
  Apk0 <- Apk(p, k0)
  Apk1 <- Apk(p, k1)
  up <- k0^(p/2 - 1) / I0 * exp(k0 * Apk0)
  down <- k1^(p/2 - 1) / I1 * exp(k1 * Apk1)
  P <- (n - 2) * ( (up / down)^( -2/(p - 1) ) - 1 )

  Pp <- numeric(B)
  for (i in 1:B) {
    ind <- Rfast2::Sample(ina, n)
    S <- rowsum(x, ind)
    Ri <- sqrt( Rfast::rowsums(S^2) )  ## the resultant length of each group
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

    I0 <- besselI( k0, p/2 - 1 )
    I1 <- besselI( k1, p/2 - 1 )
    Apk0 <- Apk(p, k0)
    Apk1 <- Apk(p, k1)
    up <- k0^(p/2 - 1) / I0 * exp(k0 * Apk0)
    down <- k1^(p/2 - 1) / I1 * exp(k1 * Apk1)
    Pp[i] <- (n - 2) * ( (up / down)^( -2/(p - 1) ) - 1 )
  }

  p.value <- ( sum(Pp > P) + 1 ) / (B + 1)
  statistic <- P  ;   names(statistic) <- "hclr test statistic"
  parameter <- "NA"     ;   names(parameter) <- "df"
  alternative <- "The 2 directional mean vector differ"
  method <- "Permutation ANOVA for 2 directional mean vectors using the high concentration log-likelihood ratio test"
  data.name <- c("data ", " groups")
  result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                  alternative = alternative, method = method, data.name = data.name )
  class(result) <- "htest"
  return(result)
}

