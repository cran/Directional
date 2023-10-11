hetcircboot <- function(u, ina, rads = TRUE, B = 999) {
  if ( !rads )  u <- u * pi/180
  x <- cbind( cos(u), sin(u) )
  mod <- Directional::hetboot(x, ina, B = B)

  statistic <- mod$statistic
  p.value <- mod$p.value
  parameter <- mod$parameter
  alternative <- "At least one circular mean differs"
  method <- "Bootstrap ANOVA for circular means using the heterogeneous approach"
  data.name <- c("data ", " groups")
  result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                  alternative = alternative, method = method, data.name = data.name )
  class(result) <- "htest"
  return(result)
}
