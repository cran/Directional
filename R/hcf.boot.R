hcf.boot <- function(x1, x2, fc = TRUE, B = 999) {

  n1 <- dim(x1)[1]    ;  n2 <- dim(x2)[1]
  x <- rbind(x1, x2)
  ina <- c( rep(1, n1), rep(2, n2) )
  ni <- c(n1, n2)
  p <- dim(x)[2]  ## dimensionality of the data
  n <- n1 + n2  ## sample size of the data
  S <- rowsum(x, ina)
  m1 <- S[1, ]   ;   m2 <- S[2, ]
  m1 <- m1 / sqrt( sum(m1^2) )
  m2 <- m2 / sqrt( sum(m2^2) )
  Ri <- sqrt( Rfast::rowsums(S^2) )  ## the resultant length of each group
  S <- Rfast::colsums(x)
  R <- sqrt( sum(S^2) )  ## the resultant length based on all the data
  ## Next we stimate the common concentration parameter kappa
  m <- S / R
  kapaa <- Directional::vmf.mle(x, fast = TRUE)$kappa
  ## kapaa is the estimated concentration parameter based on all the data
  Ft <- (n - 2) * ( sum(Ri) - R) / ( n - sum(Ri) )
  if (fc) {  ## correction is used
    if (p == 3) {
      Ft <- kapaa * (1/kapaa - 1/(5 * kapaa^3)) * Ft
    } else if (p > 3)  Ft <- kapaa * ( 1/kapaa - (p - 3)/(4 * kapaa^2) - (p - 3)/(4 * kapaa^3) ) * Ft
  }

  rot1 <- t( Directional::rotation(m1, m) )
  rot2 <- t( Directional::rotation(m2, m) )
  y1 <- x1 %*% rot1
  y2 <- x2 %*% rot2
  ftb <- numeric(B)

  for (i in 1:B) {
    b1 <- Rfast2::Sample.int(n1, n1, replace = TRUE)
    b2 <- Rfast2::Sample.int(n2, n2, replace = TRUE)
    yb <- rbind(y1[b1, ], y2[b2, ])
    S <- rowsum(yb, ina)
    Ri <- sqrt( Rfast::rowsums(S^2) )
    S <- Rfast::colsums(yb)
    R <- sqrt( sum(S^2) )
    ftb[i] <- (n - 2) * ( sum(Ri) - R) / ( n - sum(Ri) )
    if (fc) {  ## correction is used
      kapa <- Directional::vmf.mle(yb, fast = TRUE)$kappa
      if (p == 3) {
        ftb[i] <- kapa * (1/kapa - 1/(5 * kapa^3)) * ftb[i]
      } else if (p > 3)  ftb[i] <- kapa * ( 1/kapa - (p - 3)/(4 * kappa^2) - (p - 3)/(4 * kapa^3) ) * ftb[i]
    }
  }

  p.value <- ( sum(ftb > Ft) + 1 ) / (B + 1)
  statistic <- Ft  ;   names(statistic) <- "Bootstrap hcf test statistic"
  parameter <- "NA"     ;   names(parameter) <- "df"
  alternative <- "The 2 directional mean vector differ"
  method <- "Bootstrap ANOVA for 2 directional mean vectors using the high concentration test"
  data.name <- c("data ", " groups")
  result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                  alternative = alternative, method = method, data.name = data.name )
  class(result) <- "htest"
  return(result)
}
