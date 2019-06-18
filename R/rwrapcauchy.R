rwrapcauchy <- function(n, m, rho, rads = TRUE) {
  a <-  -log(rho)
  if ( !rads )  m <- m / 180 * pi
  theta <- rcauchy(n, m, a)
  theta <- theta %% (2 * pi)
  if ( !rads )  theta <- theta * pi/180
}
