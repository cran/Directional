rspcauchy <- function(n, mu, rho) {
  phi <- rho * mu
  U <- Directional::riag( n, numeric( length(mu) ) )
  com <- Rfast::eachrow(U, phi, oper = "+")
  x <- ( 1 - sum(phi^2) ) / Rfast::rowsums( com^2 ) * com
  Rfast::eachrow(x, phi, oper = "+")
}
