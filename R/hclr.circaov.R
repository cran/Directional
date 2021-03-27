hclr.circaov <- function(u, ina, rads = FALSE) {
  if ( !rads )  u <- u * pi/180
  x <- cbind( cos(u), sin(u) )
  Directional::hclr.aov(x, ina)
}

