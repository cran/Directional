################################
#### Median direction
#### Tsagris Michail 1/2016
#### mtsagris@yahoo.gr
#### References: Fisher, N. I. (1985). Spherical medians.
#### Journal of the Royal Statistical Society. Series B, 47(2): 342-348.
#### Fisher, N. I., Lewis, T., & Embleton, B. J. (1987).
#### Statistical analysis of spherical data. Cambridge university press.
################################

mediandir_2 <- function(x) {
  ## x is the directional data

  funa <- function(pa) {
    pa <- pa / sqrt( sum(pa^2) )
    mean( acos( x %*% pa ) )
  }

  bar <- nlm( funa, Rfast::colMedians(x), iterlim = 10000 )
  bar <- nlm( funa, bar$estimate, iterlim = 10000 )
  bar <- optim( bar$estimate, funa, control = list(maxit = 10000) )
  med <- bar$par
  med / sqrt( sum(med^2) )

}
