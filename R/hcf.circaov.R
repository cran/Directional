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
  Rfast2::hcf.circaov(u, ina)
}
