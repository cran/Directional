################################
#### ANOVA for cicular data (High concentration F test)
#### Tsagris Michail 1/2015
#### mtsagris@yahoo.gr
#### References: S Rao Jammalamadaka and A SenGupta (2001)
#### Topics in circular statistics, pages 125-127
#### References: Mardia Kanti V. and Jupp Peter E. (2000)
#### Directional statistics, page 135
################################
hcf.circaov <- function(u, ina, rads = FALSE) {
  ## u contains all the circular data in radians or degrees
  ## ina is an indicator variable of each sample
  ## if the data are in degrees we transform them into radians
  if ( !rads )  u <- u * pi/180
  mod <- Rfast2::hcf.circaov(u, ina)

  n <- length(u)
  ni <- tabulate(ina)
  g <- max(ina)
  statistic <- mod[1]  ;   names(statistic) <- "F-test statistic"
  p.value <- mod[2]
  parameter <- c( g - 1 , n - g )     ;   names(parameter) <- c("df1", "df2")
  alternative <- "At least one circular mean differs"
  method <- "ANOVA for circular data using the high concentration approach"
  data.name <- c("data ", " groups")
  result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                  alternative = alternative, method = method, data.name = data.name )
  class(result) <- "htest"
  return(result)
}
