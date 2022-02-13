haversine.dist <- function(x) {
  n <- dim(x)[1]
  mat <- matrix(0, nrow = n, ncol = n)
  x1 <- x[, 1]
  x2 <- x[, 2]
  for (i in 1:(n - 1) )  {
    ind <- (i + 1):n
    ind_row <- x1[ind]
    a <- sin( 0.5 * (x1[i] - ind_row) )^2 + cos(x1[i]) * cos(ind_row) * sin( 0.5 * (x2[i] - x2[ind]) )^2
    mat[i, ind] <- mat[ind, i] <- 2 * asin( sqrt(a) )
  }

  mat
}
