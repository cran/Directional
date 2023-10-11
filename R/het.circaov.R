################################
#### ANOVA for cicular data (Heterogeneous case, kappas not equal)
#### Tsagris Michail 1/2015
#### mtsagris@yahoo.gr
#### References: Mardia Kanti V. and Jupp Peter E. (2000)
#### Directional statistics, page 141-142
################################
het.circaov <- function(u, ina, rads = FALSE) {
  ## u contains all the circular data in radians or degrees
  ## ina is an indicator variable of each sample
  ## if the data are in degrees we transform them into radians
  if ( !rads )  u <- u * pi/180
  mod <- Rfast2::het.circaov(u, ina)

  n <- length(u)
  ni <- tabulate(ina)
  g <- max(ina)
  statistic <- mod[1]  ;   names(statistic) <- "chi-square test statistic"
  p.value <- mod[2]
  parameter <- g - 1      ;   names(parameter) <- "df"
  alternative <- "At least one circular mean differs"
  method <- "ANOVA for circular data using the heterogeneous approach"
  data.name <- c("data ", " groups")
  result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                  alternative = alternative, method = method, data.name = data.name )
  class(result) <- "htest"
  return(result)
}
