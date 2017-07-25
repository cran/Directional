#### Projected multivariate normal MLE
#### Presnell, Morrison and Littell (1998), JASA
################################
spml.mle <- function(x, rads = FALSE, tol = 1e-07) {
  if ( !rads )   x <- x * pi/180
  mod <- Rfast::spml.mle(x, tol = tol)
  if (!rads)   mod$mumu <- ( atan(mod$mu[2]/mod$mu[1]) + pi * I(mod$mu[1] < 0) ) %% (2 * pi)
  mod
}


