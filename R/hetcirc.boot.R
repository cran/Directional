hetcirc.boot <- function(u1, u2, rads = TRUE, B = 999) {
  u <- c(u1, u2)
  if ( !rads ) u <- u * pi/180
  x1 <- cbind( cos(u1), sin(u1) )
  x2 <- cbind( cos(u2), sin(u2) )
  Directional::het.boot(x1, x2, B = B)
}
