multispml.mle <- function(x, ina, tol = 1e-07, ell = FALSE) {
  ni <- tabulate(ina)
  ni <- ni[ni > 0]
  g <- length(ni)
  loglik <- gi <- numeric(g)
  mi <- matrix(nrow = g, ncol = 2)
  for (i in 1:g) {
    mod <- Rfast::spml.mle( x[ina == i], tol = tol )
    loglik[i] <- mod$loglik
    gi <- mod$gamma
    mi[i, ] <- mod$mu
  }
  if ( !ell )  loglik <- NULL
  list(loglik = loglik, gi = gi, mi = mi)
}













