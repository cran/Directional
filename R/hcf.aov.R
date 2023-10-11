hcf.aov <- function(x, ina, fc = TRUE) {
  ## x contains all the data
  ## ina is an indicator variable of each sample
  ina <- as.numeric(ina)
  g <- max(ina)  ## how many groups are there
  p <- dim(x)[2]
  n <- dim(x)[1]  ## dimensionality and sample size of the data
  S <- rowsum(x, ina)
  Ri <- sqrt( Rfast::rowsums(S^2) )  ## the resultant length of each group
  S <- Rfast::colsums(x)
  R <- sqrt( sum(S^2) )  ## the resultant length based on all the data
  ## Next we stimate the common concentration parameter kappa
  kappa <- Directional::vmf.mle(x, fast = TRUE)$kappa
  ## kappa is the estimated concentration parameter based on all the data
  Ft <- (n - g) * (sum(Ri) - R) /( (g - 1) * (n - sum(Ri)) )
  if (fc) {  ## correction is used
    if (p == 3) {
      Ft <- kappa * (1/kappa - 1/(5 * kappa^3)) * Ft
    } else if (p > 3)  Ft <- kappa * ( 1/kappa - (p - 3)/(4 * kappa^2) - (p - 3)/(4 * kappa^3) ) * Ft
  }

  p.value <- pf(Ft, (g - 1) * (p - 1), (n - g) * (p - 1), lower.tail = FALSE)
  statistic <- Ft  ;   names(statistic) <- "F-test statistic"
  parameter <- c( (g - 1) * (p - 1), (n - g) * (p - 1) )     ;   names(parameter) <- c("df1", "df2")
  alternative <- "At least one directional mean vector differs"
  method <- "ANOVA for directional data using the high concentration test"
  data.name <- c("data ", " groups")
  result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                  alternative = alternative, method = method, data.name = data.name )
  class(result) <- "htest"
  return(result)
}
