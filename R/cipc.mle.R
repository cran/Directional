cipc.mle <- function(x, rads = FALSE, tol = 1e-6) {

  if ( !rads )  x <- x * pi/180
  x <- cbind( cos(x), sin(x) )
  n <- dim(x)[1]
  mu <- Rfast::colmeans(x)
  g2 <- sum(mu^2)
  a <- as.vector(x %*% mu)
  com <- sqrt(g2 + 1)
  com2 <- 1 / (com - a)
  lik <- -sum( log( sqrt(g2 + 1) - a ) )
  up <- Rfast::eachrow(x, mu / com, oper = "-" )
  der <- Rfast::eachcol.apply(up, com2)
  #up1 <- Rfast::Outer( as.vector( ( diag(com, 2) - tcrossprod(mu) / com ) / com^2 ), com2, oper = "/" )
  #up1 <- matrix( Rfast::colsums(up1), ncol = 2)
  up1 <- ( diag(com, 2) - tcrossprod(mu) / com ) / com^2 * sum(com2)
  up2 <- crossprod(up * com2)
  der2 <- up2 - up1

  mu <- mu - solve(der2, der)
  g2 <- sum(mu^2)
  a <- as.vector(x %*% mu)
  com <- sqrt(g2 + 1)
  com2 <- 1 / (com - a)
  lik[2] <- sum( log( com2 ) )

  i <- 2
  while ( lik[i] - lik[i - 1] > tol ) {
    i <- i + 1
    up <- Rfast::eachrow(x, mu / com, oper = "-" )
    der <- Rfast::eachcol.apply(up, com2)
    #up1 <- Rfast::Outer( as.vector( ( diag(com, 2) - tcrossprod(mu) / com ) / com^2 ), com2, oper = "/" )
    #up1 <- matrix( Rfast::colsums(up1), ncol = 2)
    up1 <- ( diag(com, 2) - tcrossprod(mu) / com ) / com^2 * sum(com2)
    up2 <- crossprod(up * com2)
    der2 <- up2 - up1
    mu <- mu - solve(der2, der)
    a <- as.vector(x %*% mu)
    g2 <- sum(mu^2)
    com <- sqrt(g2 + 1)
    com2 <- 1 / ( com - a )
    lik[i] <- sum( log( com2 ) )
  }
  circmu <- ( atan(mu[2]/mu[1]) + pi * I(mu[1] < 0) ) %% (2 * pi)
  gamma <- sqrt( sum(mu^2) )
  list(mu = mu, circmu = circmu, gamma = gamma, loglik = lik[i] - n * log(2 * pi) )
}







# cipc.mle <- function(x, rads = FALSE) {
#
#   lik <- function(mu, x) {
#     g2 <- sum(mu^2)
#     a <- as.vector(x %*% mu)
#     sum( log( sqrt(g2 + 1) - a ) )
#   }
#
#   if ( !rads )  x <- x * pi/180
#   x <- cbind( cos(x), sin(x) )
#   n <- dim(x)[1]
#
#   mod <- optim( rnorm(2), lik, x = x, control = list(maxit = 5000) )
#   mod <- optim( mod$par, lik, x = x, control = list(maxit = 5000) )
#   mu <- mod$par
#   circmu <- ( atan(mu[2]/mu[1]) + pi * I(mu[1] < 0) ) %% (2 * pi)
#   gama <- sqrt( sum(mu^2) )
#   list(mu = mu, circmu = circmu, gamma = gama, loglik = -mod$value - n * log(2 * pi) )
# }
