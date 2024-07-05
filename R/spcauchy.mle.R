spcauchy.mle <- function(x, tol = 1e-6) {

  dm <- dim(x)
  n <- dm[1]  ;  d <- dm[2] - 1
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
  up1 <- ( diag(com, d + 1) - tcrossprod(mu) / com ) / com^2 * sum(com2)
  up2 <- crossprod(up * com2)
  der2 <- up2 - up1

  mu <- mu - solve(der2, der)
  g2 <- sum(mu^2)
  a <- as.vector(x %*% mu)
  com <- sqrt(g2 + 1)
  com2 <- 1 / (com - a)
  lik[2] <- d * sum( log( com2 ) )

  i <- 2
  while ( lik[i] - lik[i - 1] > tol ) {
    i <- i + 1
    up <- Rfast::eachrow(x, mu / com, oper = "-" )
    der <- Rfast::eachcol.apply(up, com2)
    #up1 <- Rfast::Outer( as.vector( ( diag(com, 2) - tcrossprod(mu) / com ) / com^2 ), com2, oper = "/" )
    #up1 <- matrix( Rfast::colsums(up1), ncol = 2)
    up1 <- ( diag(com, d + 1) - tcrossprod(mu) / com ) / com^2 * sum(com2)
    up2 <- crossprod(up * com2)
    der2 <- up2 - up1
    mu <- mu - solve(der2, der)
    a <- as.vector(x %*% mu)
    g2 <- sum(mu^2)
    com <- sqrt(g2 + 1)
    com2 <- 1 / ( com - a )
    lik[i] <- d * sum( log( com2 ) )
  }
  gamma <- sqrt( sum(mu^2) )
  list( mesos = mu, mu = mu / gamma, gamma = gamma, rho = (com - 1) / gamma, loglik = lik[i] +
        n * lgamma( 0.5 * (d + 1) ) - 0.5 * n * (d + 1) * log(pi) - n * log(2) )
}





