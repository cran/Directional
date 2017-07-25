################################
#### Median direction
#### Tsagris Michail 1/2016
#### mtsagris@yahoo.gr
#### References: Fisher, N. I. (1985). Spherical medians.
#### Journal of the Royal Statistical Society. Series B, 47(2): 342-348.
#### Fisher, N. I., Lewis, T., & Embleton, B. J. (1987).
#### Statistical analysis of spherical data. Cambridge university press.
#### Cabrera, J., & Watson, G. S. (1990). On a spherical median related distribution.
#### Communications in Statistics-Theory and Methods, 19(6): 1973-1986.
################################
mediandir <- function(x) {
  Rfast::mediandir(x)
}
