lrcirc.boot <- function(u1, u2, rads = TRUE, B = 999) {
  if ( !rads )  {
    u1 <- u1 * pi/180
    u2 <- u2 * pi/180
  }
  u <- c(u1, u2)
  n1 <- length(u1)   ;  n2 <- length(u2)
  n <- n1 + n2
  ina <- c( rep(1, n1), rep(2, n2) )
  ni <- c(n1, n2)
  x <- cbind(cos(u), sin(u))
  rsi <- rowsum(x, ina)
  Ri <- sqrt(Rfast::rowsums(rsi^2))
  ni <- tabulate(ina)
  mi <- rsi/ni
  mi <- mi/sqrt(Rfast::rowsums(mi^2))
  m <- Rfast::colmeans(x)
  rs <- n * m
  m <- m/sqrt(sum(m^2))
  m <- matrix(rep(m, 2), nrow = 2, byrow = TRUE)

  mu <- atan(rs[2]/rs[1]) + pi * (rs[1] < 0)
  con <- sum( cos(u - mu) )
  R <- sqrt( sum(rs^2) )
  k1 <- (1.28 - 0.53 * R^2/n^2) * tan(0.5 * pi * R/n)
  if (k1 < 710) {
    der <- con - n * besselI(k1, 1, expon.scaled = TRUE) / besselI(k1, 0, expon.scaled = TRUE)
    a <- besselI(k1, 0)^2/2 + besselI(k1, 2) * besselI(k1, 0)/2 - besselI(k1, 1)^2
    der2 <- n * a/besselI(k1, 0)^2
    k2 <- k1 + der/der2
    while ( abs(k1 - k2) > 1e-08 ) {
      k1 <- k2
      der <- con - n * besselI(k1, 1, expon.scaled = TRUE) / besselI(k1, 0, expon.scaled = TRUE)
      a <- besselI(k1, 0)^2/2 + besselI(k1, 2) * besselI(k1, 0)/2 - besselI(k1, 1)^2
      der2 <- n * a/besselI(k1, 0)^2
      k2 <- k1 + der/der2
    }
  } else k2 <- k1
  kapaa <- k2
  w <- kapaa * sum( Ri * Rfast::rowsums( (mi - m)^2 ) )

  m1 <- mi[1, ]
  m2 <- mi[2, ]
  rot1 <- t( Directional::rotation(m1, m[1, ]) )
  rot2 <- t( Directional::rotation(m2, m[1, ]) )
  y1 <- x[1:n1, ] %*% rot1
  y2 <- x[-c(1:n1), ] %*% rot2
  wb <- numeric(B)
  for (i in 1:B) {
    b1 <- sample(n1, n1, replace = TRUE)
    b2 <- sample(n2, n2, replace = TRUE)
    yb <- rbind(y1[b1, ], y2[b2, ])
    rsi <- rowsum(yb, ina)
    Ri <- sqrt(Rfast::rowsums(rsi^2))
    mi <- rsi/ni
    mi <- mi/sqrt(Rfast::rowsums(mi^2))
    m <- Rfast::colmeans(yb)
    rs <- n * m
    m <- m/sqrt(sum(m^2))
    m <- matrix(rep(m, 2), nrow = 2, byrow = TRUE)

    mu <- atan(rs[2]/rs[1]) + pi * (rs[1] < 0)
    ub <- atan(yb[, 2]/yb[, 1]) + pi * (yb[, 1] < 0)
    con <- sum( cos(ub - mu) )
    R <- sqrt( sum(rs^2) )
    k1 <- (1.28 - 0.53 * R^2/n^2) * tan(0.5 * pi * R/n)
    if (k1 < 710) {
      der <- con - n * besselI(k1, 1, expon.scaled = TRUE) / besselI(k1, 0, expon.scaled = TRUE)
      a <- besselI(k1, 0)^2/2 + besselI(k1, 2) * besselI(k1, 0)/2 - besselI(k1, 1)^2
      der2 <- n * a/besselI(k1, 0)^2
      k2 <- k1 + der/der2
      while ( abs(k1 - k2) > 1e-08 ) {
        k1 <- k2
        der <- con - n * besselI(k1, 1, expon.scaled = TRUE) / besselI(k1, 0, expon.scaled = TRUE)
        a <- besselI(k1, 0)^2/2 + besselI(k1, 2) * besselI(k1, 0)/2 - besselI(k1, 1)^2
        der2 <- n * a/besselI(k1, 0)^2
        k2 <- k1 + der/der2
      }
    } else k2 <- k1
    kapa <- k2
    wb[i] <- kapa * sum( Ri * Rfast::rowsums( (mi - m)^2 ) )
  }  ## end for (i in 1:B)

  pvalue <- ( sum(wb > w) + 1 ) / (B + 1)
  res <- c(w, pvalue, kapaa)
  names(res) <- c("test", "p-value", "kappa")
  res
}
