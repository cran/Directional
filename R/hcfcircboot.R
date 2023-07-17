hcfcircboot <- function(u, ina, rads = TRUE, B = 999) {
  if ( !rads )  u <- u * pi/180
  x <- cbind( cos(u), sin(u) )
  Directional::hcfboot(x, ina, B = B)
}
