purka.mle <- function(x, tol = 1e-07) {
  mod <- Rfast2::purka.mle(x, tol = tol)
  if ( length( mod$theta ) == 2 ) {
    m <- mod$theta
    circtheta <- ( atan(m[2] / m[1]) + pi * I(m[1] < 0) ) %%(2 * pi)
    res <- list(theta = mod$theta, circtheta = circtheta, alpha = mod$alpha, loglik = mod$loglik, alpha.sd = mod$alpha.sd)
  } else  res <- mod
  res
}





