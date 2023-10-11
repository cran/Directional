het.boot <- function(x1, x2, B = 999) {

  n1 <- dim(x1)[1]    ;  n2 <- dim(x2)[1]
  x <- rbind(x1, x2)
  ina <- c( rep(1, n1), rep(2, n2) )
  ni <- c(n1, n2)
  p <- dim(x)[2]
  n <- n1 + n2
  kapa <- numeric(2)
  mi <- rowsum(x, ina) / ni
  Ri <- sqrt( Rfast::rowsums(mi^2) )
  m <- mi/Ri
  m1 <- m[1, ]    ;   m2 <- m[2, ]
  S <- Rfast::colsums(x)
  R <- sqrt( sum(S^2) )
  m <- S/R

  kapa[1] <- Directional::vmf.mle( x1, fast = TRUE )$kappa
  kapa[2] <- Directional::vmf.mle( x2, fast = TRUE )$kappa
  tw <- Rfast::colsums(kapa * ni * mi)
  Tt <- sum( kapa * ni * sqrt( Rfast::rowsums(mi^2) ) ) - sqrt( sum(tw^2) )

  tb <- numeric(B)
  rot1 <- t( Directional::rotation(m1, m) )
  rot2 <- t( Directional::rotation(m2, m) )
  y1 <- x1 %*% rot1
  y2 <- x2 %*% rot2

  for (i in 1:B) {

    b1 <- Rfast2::Sample.int(n1, n1, replace = TRUE)
    b2 <- Rfast2::Sample.int(n2, n2, replace = TRUE)
    yb1 <- y1[b1, ]   ;   yb2 <- y2[b2, ]
    yb <- rbind(yb1, yb2)
    mi <- rowsum(yb, ina) / ni
    kapa[1] <- Directional::vmf.mle(yb1, fast = TRUE )$kappa
    kapa[2] <- Directional::vmf.mle(yb2, fast = TRUE )$kappa
    tw <- Rfast::colsums(kapa * ni * mi)
    tb[i] <- sum( kapa * ni * sqrt( Rfast::rowsums(mi^2) ) ) - sqrt( sum(tw^2) )
  }

  p.value <- ( sum(tb > Tt) + 1 ) / (B + 1)
  statistic <- 2 * Tt  ;   names(statistic) <- "Bootstrap het test statistic"
  parameter <- "NA"     ;   names(parameter) <- "df"
  alternative <- "The 2 directional mean vectors differ"
  method <- "Bootstrap ANOVA for 2 directional mean vectors using the heterogeneous approach"
  data.name <- c("data ", " groups")
  result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                  alternative = alternative, method = method, data.name = data.name )
  class(result) <- "htest"
  return(result)
}
