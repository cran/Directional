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
    a <- x %*% pa
    a[ abs(a) > 1 ] <- 1
    mean( acos( a ) )
  }
  options(warn = -1)
  bar <- nlm( funa, Rfast::colMedians(x), iterlim = 10000 )
  ini <- bar$estimate/sqrt( sum(bar$estimate^2) )
  bar <- optim( ini, funa, control = list(maxit = 10000) )
  med <- bar$par
  med / sqrt( sum(med^2) )
}
