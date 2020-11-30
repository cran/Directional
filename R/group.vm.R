group.vm <- function(group, fi, rads = FALSE) {
  group <- as.matrix(group)
  if ( !rads )   group <- group / 180 * pi
  u <- Rfast::rowmeans(group) ## mid points of the classes
  u <- rep(u, fi)
  mod <- Directional::circ.summary(u, rads = TRUE, plot = FALSE)
  h <- mean( apply(group, 1, diff) ) ## mean range of the classes
  ah <-  0.5 * h / sin( 0.5 * h )
  mod$MRL <- ah * mod$MRL ## grouped data correction of the MRL
  mod
}
