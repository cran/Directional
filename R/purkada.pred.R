purkada.pred <- function(ynew, y, ina) {
  ## xnew is the new observation(s)
  ## x is the data set
  ## ina is the group indicator variable
  ynew <- as.matrix(ynew)
  if ( ncol(ynew) == 1 )   ynew <- t(ynew)
  ina <- as.numeric(ina)
  g <- max(ina)
  mat <- matrix(0, dim(ynew)[1], g)
  for (j in 1:g) {
    mod <- Directional::purka.mle( y[ina == j, ] )
    mat[, j] <- Directional::dpurka(ynew, mod$theta, mod$alpha, logden = TRUE )
  }
  Rfast::rowMaxs(mat)
}
