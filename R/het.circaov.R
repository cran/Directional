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
  Rfast2::het.circaov(u, ina)
}
