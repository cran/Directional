################################
#### ANOVA for cicular data (Embedding approach)
#### Tsagris Michail 1/2015
#### mtsagris@yahoo.gr
#### References: Mardia Kanti V. and Jupp Peter E. (2000)
#### Directional statistics, page 138-139
################################
embed.circaov <- function(u, ina, rads = FALSE) {
  if ( !rads )  u <- u * pi/180
  Rfast2::embed.circaov(u, ina)
}
