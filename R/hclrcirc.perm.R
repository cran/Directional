hclrcirc.perm <- function(u1, u2, rads = TRUE, B = 999) {

  if ( !rads )  {
    u1 <- u1 * pi/180
    u2 <- u2 * pi/180
  }
  x1 <- cbind( cos(u1), sin(u1) )
  x2 <- cbind( cos(u2), sin(u2) )
  Directional::hclr.perm(x1, x2)
}
