hetboot <- function(x, ina, B = 999) {

  ina <- as.numeric(ina)
  x <- x[order(ina), ]
  ina <- sort(ina)
  ni <- tabulate(ina)
  k <- length(ni)
  dm <- dim(x)
  p <- dm[2]  ## dimensionality of the data
  n <- dm[1]  ## sample size of the data

  kapa <- numeric(k)
  mi <- rowsum(x, ina) / ni
  Ri <- sqrt( Rfast::rowsums(mi^2) )
  S <- Rfast::colsums(x)
  R <- sqrt( sum(S^2) )
  m <- S/R

  for (j in 1:k)  kapa[j] <- Directional::vmf.mle( x[ina == j, ], fast = TRUE )$kappa
  tw <- Rfast::colsums(kapa * ni * mi)
  Tt <- sum( kapa * ni * sqrt( Rfast::rowsums(mi^2) ) ) - sqrt( sum(tw^2) )

  mi <- mi/Ri
  y <- list()
  for (j in 1:k) {
    rot <- t( Directional::rotation(mi[j, ], m) )
    y[[ j ]] <- x[ina == j, ] %*% rot
  }
  tb <- numeric(B)

  for (i in 1:B) {

    yb <- NULL
    for (j in 1:k) {
      b <- Rfast2::Sample.int(ni[j], ni[j], replace = TRUE)
      yb <- rbind( yb, y[[ j ]][b, ] )
      kapa[j] <- Directional::vmf.mle( y[[ j ]][b, ], fast = TRUE )$kappa
    }
    mi <- rowsum(yb, ina) / ni
    tw <- Rfast::colsums(kapa * ni * mi)
    tb[i] <- sum( kapa * ni * sqrt( Rfast::rowsums(mi^2) ) ) - sqrt( sum(tw^2) )
  }

  p.value <- ( sum(tb > Tt) + 1 ) / (B + 1)
  statistic <- 2 * Tt  ;   names(statistic) <- "Bootstrap het test statistic"
  parameter <- "NA"     ;   names(parameter) <- "df"
  alternative <- "At least one directional mean vector differs"
  method <- "Bootstrap ANOVA for directional data using the heterogeneous approach"
  data.name <- c("data ", " groups")
  result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                  alternative = alternative, method = method, data.name = data.name )
  class(result) <- "htest"
  return(result)
}

