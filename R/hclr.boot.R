hclr.boot <- function(x1, x2, B = 999) {

  n1 <- dim(x1)[1]    ;  n2 <- dim(x2)[1]
  x <- rbind(x1, x2)
  ina <- c( rep(1, n1), rep(2, n2) )
  ni <- c(n1, n2)
  p <- dim(x)[2]
  n <- n1 + n2
  S <- rowsum(x, ina)
  Ri <- sqrt( Rfast::rowsums(S^2) )
  m <- S/Ri
  m1 <- m[1, ]    ;   m2 <- m[2, ]
  S <- Rfast::colsums(S)
  R <- sqrt( sum(S^2) )
  m <- S/R
  Apk <- function(p, k)  besselI(k, p/2, expon.scaled = TRUE) / besselI(k, p/2 - 1, expon.scaled = TRUE)

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

  Pb <- numeric(B)
  rot1 <- t( Directional::rotation(m1, m) )
  rot2 <- t( Directional::rotation(m2, m) )
  y1 <- x1 %*% rot1
  y2 <- x2 %*% rot2

  for (i in 1:B) {

    b1 <- Rfast2::Sample.int(n1, n1, replace = TRUE)
    b2 <- Rfast2::Sample.int(n2, n2, replace = TRUE)
    yb <- rbind(y1[b1, ], y2[b2, ])
    S <- rowsum(yb, ina)
    Ri <- sqrt( Rfast::rowsums(S^2) )
    S <- Rfast::colsums(yb)
    R <- sqrt( sum(S^2) )

    Rk <- R/n
    k1 <- Rk * (p - Rk^2)/(1 - Rk^2)
    k2 <- k1 - (Apk(p, k1) - Rk) / ( 1 - Apk(p, k1)^2 - (p - 1)/k1 * Apk(p, k1) )
    while (abs(k2 - k1) > 1e-07) {
      k1 <- k2
      k2 <- k1 - (Apk(p, k1) - Rk) / (1 - Apk(p, k1)^2 - (p - 1)/k1 * Apk(p, k1) )
    }
    k0 <- k2  ## concentration parameter under H0

    Rk <- sum(Ri)/n
    k1 <- Rk * (p - Rk^2)/(1 - Rk^2)
    k2 <- k1 - (Apk(p, k1) - Rk) / ( 1 - Apk(p, k1)^2 - (p - 1)/k1 * Apk(p, k1) )

    while (abs(k2 - k1) > 1e-07) {
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
    Pb[i] <- (n - 2) * ( (up / down)^( -2/(p - 1) ) - 1 )

  }

  p.value <- ( sum(Pb > P) + 1 ) / (B + 1)
  statistic <- P  ;   names(statistic) <- "Bootstrap hclr test statistic"
  parameter <- "NA"     ;   names(parameter) <- "df"
  alternative <- "The 2 directional mean vector differ"
  method <- "Bootstrap ANOVA for 2 directional mean vectors using the high concentration log-likelihood ratio test"
  data.name <- c("data ", " groups")
  result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                  alternative = alternative, method = method, data.name = data.name )
  class(result) <- "htest"
  return(result)
}

