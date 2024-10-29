dmixspcauchy <- function(y, probs, mu, rho, logden = FALSE) {
  g <- length(rho)
  den <- matrix(nrow = dim(y)[1], ncol = g)
  for (j in 1:g)  den[, j] <- Directional::dspcauchy(y, mu[j, ], rho[j])
  den <- den %*% probs
  if ( logden )  den <- log(den)
  den
}
