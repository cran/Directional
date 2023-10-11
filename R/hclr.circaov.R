hclr.circaov <- function(u, ina, rads = FALSE) {
  if ( !rads )  u <- u * pi/180
  x <- cbind( cos(u), sin(u) )
  mod <- Directional::hclr.aov(x, ina)

  statistic <- mod$statistic
  p.value <- mod$p.value
  parameter <- mod$parameter
  alternative <- "At least one circular mean differs"
  method <- "ANOVA for circular data using the high concentration log-likelihood ratio test"
  data.name <- c("data ", " groups")
  result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                  alternative = alternative, method = method, data.name = data.name )
  class(result) <- "htest"
  return(result)
}

