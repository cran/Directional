hcf.perm <- function(x1, x2, B = 999) {

  n1 <- dim(x1)[1]    ;  n2 <- dim(x2)[1]
  x <- rbind(x1, x2)
  ina <- c( rep(1, n1), rep(2, n2) )
  ni <- c(n1, n2)
  p <- dim(x)[2]  ## dimensionality of the data
  n <- n1 + n2  ## sample size of the data
  S <- rowsum(x, ina)
  Ri <- sqrt( Rfast::rowsums(S^2) )  ## the resultant length of each group
  S <- Rfast::colsums(x)
  R <- sqrt( sum(S^2) )  ## the resultant length based on all the data

  Ft <- (n - 2) * (sum(Ri) - R) / ( n - sum(Ri) )

  pft <- numeric(B)
  for (i in 1:B) {
    ind <- Rfast2::Sample(ina, n)
    S <- rowsum(x, ind)
    Ri <- sqrt( Rfast::rowsums(S^2) )
    pft[i] <- (n - 2) * (sum(Ri) - R) / ( n - sum(Ri) )
  }

  p.value <- ( sum(pft > Ft) + 1 ) / (B + 1)
  statistic <- Ft  ;   names(statistic) <- "hcf test statistic"
  parameter <- "NA"     ;   names(parameter) <- "df"
  alternative <- "The 2 directional mean vector differ"
  method <- "Permutation ANOVA for 2 directional mean vectors using the high concentration test"
  data.name <- c("data ", " groups")
  result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                  alternative = alternative, method = method, data.name = data.name )
  class(result) <- "htest"
  return(result)
}
