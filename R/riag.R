riag <- function(n, mu) {

  y <- Rfast::matrnorm( n, length(mu) )
  y <- Rfast::eachrow(y, mu, oper = "+")
  y / sqrt( Rfast::rowsums(y^2) )

}
