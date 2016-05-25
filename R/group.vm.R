group.vm <- function(group, fi, rads = FALSE) {
  
  group <- as.matrix(group) 
  if ( rads == FALSE )   group <- group / 180 * pi 
  u <- rowMeans(group) ## mid points of the classes
  u <- rep(u, fi)

  mod <- circ.summary(u, rads = TRUE, plot = FALSE)
  h <- mean( apply(group, 1, diff) ) ## mean range of the classes
  ah <- ( h / 2 ) / ( sin( h / 2 ) )
  mod$MRL <- ah * mod$MRL ## grouped data correction of the MRL
  mod 

} 