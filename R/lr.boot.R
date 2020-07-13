lr.boot <- function(x1, x2, B = 999) {

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
  S <- Rfast::colsums(x)
  R <- sqrt( sum(S^2) )
  m <- S/R
  Apk <- function(p, k)  besselI(k, p/2, expon.scaled = TRUE) / besselI(k, p/2 - 1, expon.scaled = TRUE)

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

  apk0 <- (1 - p/2) * log(k0/2) + lgamma(p/2) + log( besselI( k0, p/2 - 1, expon.scaled = TRUE ) ) + k0
  apk1 <- (1 - p/2) * log(k1/2) + lgamma(p/2) + log( besselI( k1, p/2 - 1, expon.scaled = TRUE ) ) + k1
  w <- 2 * (k1 * sum(Ri) - k0 * R - n * apk1 + n * apk0)

  wb <- numeric(B)
  rot1 <- t( Directional::rotation(m1, m) )
  rot2 <- t( Directional::rotation(m2, m) )
  y1 <- x1 %*% rot1
  y2 <- x2 %*% rot2

  for (i in 1:B) {

    b1 <- sample(n1, n1, replace = TRUE)
    b2 <- sample(n2, n2, replace = TRUE)
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

    apk0 <- (1 - p/2) * log(k0/2) + lgamma(p/2) + log( besselI( k0, p/2 - 1, expon.scaled = TRUE ) ) + k0
    apk1 <- (1 - p/2) * log(k1/2) + lgamma(p/2) + log( besselI( k1, p/2 - 1, expon.scaled = TRUE ) ) + k1
    wb[i] <- 2 * (k1 * sum(Ri) - k0 * R - n * apk1 + n * apk0)

  }

  pvalue <- ( sum(wb > w) + 1 ) / (B + 1)
  res <- c(w, pvalue)
  names(res) <- c('w', 'p-value')
  res
}

