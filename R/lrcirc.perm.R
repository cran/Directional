lrcirc.perm <- function(u1, u2, rads = TRUE, B = 999) {

  u <- c(u1, u2)
  if ( !rads )  u <- u * pi/180
  ina <- c( rep(1, length(u1) ), rep(2, length(u2) ) )
  ni <- tabulate(ina)
  x <- cbind( cos(u), sin(u) )
  n <- dim(x)[1]

  rsi <- rowsum(x, ina)
  Ri <- sqrt( Rfast::rowsums( rsi^2) )
  mi <- rsi/ni
  mi <- mi/sqrt( Rfast::rowsums(mi^2) )
  m <- Rfast::colmeans(x)
  m <- m/sqrt( sum(m^2) )
  m <- matrix( rep(m, 2), nrow = 2, byrow = TRUE )

  rs <- Rfast::colsums(rsi)
  mu <- atan( rs[2]/rs[1] ) + pi * (rs[1] < 0)
  con <- sum( cos(u - mu) )
  R <- sqrt( sum(rs^2) )
  k1 <- (1.28 - 0.53 * R^2/n^2) * tan(0.5 * pi * R/n)
  if (k1 < 710) {
    der <- con - n * besselI(k1, 1, expon.scaled = TRUE)/besselI(k1, 0, expon.scaled = TRUE)
    a <- besselI(k1, 0)^2/2 + besselI(k1, 2) * besselI(k1, 0)/2 - besselI(k1, 1)^2
    der2 <- n * a/besselI(k1, 0)^2
    k2 <- k1 + der/der2
    while (abs(k1 - k2) > 1e-07) {
      k1 <- k2
      der <- con - n * besselI(k1, 1, expon.scaled = TRUE)/besselI(k1, 0, expon.scaled = TRUE)
      a <- besselI(k1, 0)^2/2 + besselI(k1, 2) * besselI(k1, 0)/2 - besselI(k1, 1)^2
      der2 <- n * a/besselI(k1, 0)^2
      k2 <- k1 + der/der2
    }
  } else k2 <- k1
  kapa <- k2
  w <- kapa * sum( Ri * Rfast::rowsums( (mi - m)^2 ) )

  wp <- numeric(B)
  for (i in 1:B) {
    ind <- sample(ina, n)
    rsi <- rowsum(x, ind)
    Ri <- sqrt( Rfast::rowsums(rsi^2) )
    mi <- rsi/ni
    mi <- mi/sqrt( Rfast::rowsums(mi^2) )

    rs <- Rfast::colsums(rsi)
    mu <- atan( rs[2]/rs[1] ) + pi * (rs[1] < 0)
    con <- sum( cos(u - mu) )
    R <- sqrt( sum(rs^2) )
    k1 <- (1.28 - 0.53 * R^2/n^2) * tan(0.5 * pi * R/n)
    if (k1 < 710) {
      der <- con - n * besselI(k1, 1, expon.scaled = TRUE)/besselI(k1, 0, expon.scaled = TRUE)
      a <- besselI(k1, 0)^2/2 + besselI(k1, 2) * besselI(k1, 0)/2 - besselI(k1, 1)^2
      der2 <- n * a/besselI(k1, 0)^2
      k2 <- k1 + der/der2
      while (abs(k1 - k2) > 1e-07) {
        k1 <- k2
        der <- con - n * besselI(k1, 1, expon.scaled = TRUE)/besselI(k1, 0, expon.scaled = TRUE)
        a <- besselI(k1, 0)^2/2 + besselI(k1, 2) * besselI(k1, 0)/2 - besselI(k1, 1)^2
        der2 <- n * a/besselI(k1, 0)^2
        k2 <- k1 + der/der2
      }
    } else k2 <- k1
    kapa <- k2
    wp[i] <- kapa * sum( Ri * Rfast::rowsums( (mi - m)^2 ) )
  }

  pvalue <- ( sum(wp > w) + 1 ) / (B + 1)
  res <- c(kapa * w, pvalue, kapa)
  names(res) <- c("test", "p-value", "kappa")
  res
}
