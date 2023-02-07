dmixvmf <- function(y, probs, mu, k, logden = FALSE) {
  g <- length(k)
  den <- matrix(nrow = dim(y)[1], ncol = g)
  for (j in 1:g)  den[, j] <- Directional::dvmf(y, mu[j, ], k[j])
  den <- den %*% probs
  if ( logden )  den <- log(den)
  den
}
