riag <- function(n, mu) {

  y <- rangen::Rnorm.mat( n, length(mu) )
  y <- Rfast::eachrow(y, mu, oper = "+")
  y / sqrt( Rfast::rowsums(y^2) )

}
