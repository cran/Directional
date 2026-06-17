rcircexp <- function(n, lambda, rads = TRUE) {
  theta <- rangen::Rexp(n, lambda)
  theta <- theta %% (2 * pi)
  if ( !rads )  theta <- theta * pi/180
  theta
}
