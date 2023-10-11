embedcirc.boot <- function(u1, u2, rads = TRUE, B = 999) {
  if ( !rads )  {
    u1 <- u1 * pi/180
    u2 <- u2 * pi/180
  }
  x1 <- cbind( cos(u1), sin(u1) )
  x2 <- cbind( cos(u2), sin(u2) )
  mod <- Directional::embed.boot(x1, x2, B = B)

  statistic <- mod$statistic
  p.value <- mod$p.value
  parameter <- mod$parameter
  alternative <- "The 2 circular means differ"
  method <- "Bootstrap ANOVA for 2 circular means using the embedding approach"
  data.name <- c("data ", " groups")
  result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                  alternative = alternative, method = method, data.name = data.name )
  class(result) <- "htest"
  return(result)
 }
