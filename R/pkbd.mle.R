pkbd.mle <- function(x, tol = 1e-6) {

  dm <- dim(x)
  n <- dm[1]  ;  d <- dm[2] - 1
  mu <- colmeans(x)
  g2 <- sum(mu^2)
  a <- as.vector(x %*% mu)
  com <- sqrt(g2 + 1)
  com2 <- 1 / (com - a)
  com3 <- 1 / (com - 1)
  lik <- 0.5 * (d + 1) * sum( log( com2 ) ) - 0.5 * n * (d - 1) * ( log(com - 1) - log(g2) )
  up <- Rfast::eachrow(x, mu / com, oper = "-" )
  der <-  0.5 * (d + 1) * Rfast::eachcol.apply(up, com2) -
    0.5 * n * (d - 1) * ( mu / com / (com - 1) - 2 * mu / g2 )
  tcmu <- tcrossprod(mu)
  upa1 <- ( diag(com, d + 1) - tcmu / com ) / com^2 * sum(com2)
  upa2 <- crossprod(up * com2)
  #der2 <- upa2 - upa1
  upb1 <- ( diag(com, d + 1) - tcmu / com ) / com^2 * com3
  upb2 <- tcrossprod(mu) / com^2 * com3^2
  last <- ( 2 * diag(g2, d + 1) - 4 * tcmu ) / g2^2
  der2 <- 0.5 * (d + 1) * ( upa2 - upa1 ) - n * 0.5 * (d - 1) * ( upb1 - upb2 - last)
  mu <- mu - solve(der2, der)

  g2 <- sum(mu^2)
  a <- as.vector(x %*% mu)
  com <- sqrt(g2 + 1)
  com2 <- 1 / (com - a)
  lik[2] <- 0.5 * (d + 1) * sum( log( com2 ) ) - 0.5 * n * (d - 1) * ( log(com - 1) - log(g2) )

  i <- 2
  while ( abs( lik[i] - lik[i - 1] ) > tol ) {
    i <- i + 1
    up <- Rfast::eachrow(x, mu / com, oper = "-" )
    der <-  0.5 * (d + 1) * Rfast::eachcol.apply(up, com2) -
      0.5 * n * (d - 1) * ( mu / com / (com - 1) - 2 * mu / g2 )
    tcmu <- tcrossprod(mu)
    upa1 <- ( diag(com, d + 1) - tcmu / com ) / com^2 * sum(com2)
    upa2 <- crossprod(up * com2)
    #der2 <- upa2 - upa1
    upb1 <- ( diag(com, d + 1) - tcmu / com ) / com^2 * com3
    upb2 <- tcrossprod(mu) / com^2 * com3^2
    last <- ( 2 * diag(g2, d + 1) - 4 * tcmu ) / g2^2
    der2 <- 0.5 * (d + 1) * ( upa2 - upa1 ) - n * 0.5 * (d - 1) * ( upb1 - upb2 - last)
    mu <- mu - solve(der2, der)

    a <- as.vector(x %*% mu)
    g2 <- sum(mu^2)
    com <- sqrt(g2 + 1)
    com2 <- 1 / ( com - a )
    com3 <- 1 / (com - 1)
    lik[i] <-  0.5 * (d + 1) * sum( log( com2 ) ) - 0.5 * n * (d - 1) * ( log(com - 1) - log(g2) )
  }
  g <- sqrt(g2)

  list( mesos = mu, mu = mu / g, gamma = g, rho = (com - 1) / g, loglik = lik[i] - 0.5 * n * (d - 1) * log(2) +
          n * lgamma(0.5 * (d + 1) ) - n * 0.5 * (d + 1) * log(pi) - n * log(2) )
}
