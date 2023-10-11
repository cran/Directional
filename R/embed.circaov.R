################################
#### ANOVA for cicular data (Embedding approach)
#### Tsagris Michail 1/2015
#### mtsagris@yahoo.gr
#### References: Mardia Kanti V. and Jupp Peter E. (2000)
#### Directional statistics, page 138-139
################################
embed.circaov <- function(u, ina, rads = FALSE) {
  if ( !rads )  u <- u * pi/180
  mod <- Rfast2::embed.circaov(u, ina)

  n <- length(u)
  ni <- tabulate(ina)
  g <- max(ina)
  statistic <- mod[1]  ;   names(statistic) <- "F-test statistic"
  p.value <- mod[2]
  parameter <- c( g - 1 , n - g )     ;   names(parameter) <- c("df1", "df2")
  alternative <- "At least one circular mean differs"
  method <- "ANOVA for circular data using the embedding approach"
  data.name <- c("data ", " groups")
  result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                  alternative = alternative, method = method, data.name = data.name )
  class(result) <- "htest"
  return(result)
}
