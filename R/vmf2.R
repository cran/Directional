vmf2 <- function(y1, y2, tol = 1e-6) {

  dm <- dim(y1)
  n1 <- dm[1]  ;  p <- dm[2] 
  n2 <- dim(y2)[1]
  n <- n1 + n2
  sy1 <- Rfast::colsums(y1)  ;  sy2 <- Rfast::colsums(y2)
  mu <- (sy1 + sy2)
  mu <- mu / sqrt( sum(mu^2 ) )  
  
  vmfa <- function(k1, k2, R1, R2, n1, n2, p) {
    n1 * ( (p/2 - 1) * log(k1) - log(besselI(k1, p/2 - 1, expon.scaled = TRUE)) - k1 ) + k1 * R1 + 
    n2 * ( (p/2 - 1) * log(k2) - log(besselI(k2, p/2 - 1, expon.scaled = TRUE)) - k2 ) + k2 * R2       
  }

  vmfb <- function(k2, k1, R1, R2, n1, n2, p) {
    n1 * ( (p/2 - 1) * log(k1) - log(besselI(k1, p/2 - 1, expon.scaled = TRUE)) - k1 ) + k1 * R1 + 
    n2 * ( (p/2 - 1) * log(k2) - log(besselI(k2, p/2 - 1, expon.scaled = TRUE)) - k2 ) + k2 * R2 
  }

  R1 <- sum( as.vector(y1 %*% mu) )
  R2 <- sum( as.vector(y2 %*% mu) )
	
  mod1 <- optimize( vmfa, c(0, 7000), k2 = 1, R1 = R1, R2 = R2, n1 = n1, n2 = n2, p = p,
                    maximum = TRUE, tol = 1e-6 )
  mod2 <- optimize( vmfb, c(0, 7000), k1 = mod1$maximum, R1 = R1, R2 = R2, n1 = n1, n2 = n2, p = p,
                    maximum = TRUE, tol = 1e-6 )
  lik1 <- mod2$objective
  k1 <- mod1$maximum   ;   k2 <- mod2$maximum
  mu <- k1 * sy1 + k2 * sy2
  mu <- mu / sqrt( sum(mu^2) )
  R1 <- sum( as.vector(y1 %*% mu) )
  R2 <- sum( as.vector(y2 %*% mu) )

  mod1 <- optimize( vmfa, c(0, 7000), k2 = k2, R1 = R1, R2 = R2, n1 = n1, n2 = n2, p = p,
                    maximum = TRUE, tol = 1e-6 )
  mod2 <- optimize( vmfb, c(0, 7000), k1 = mod1$maximum, R1 = R1, R2 = R2, n1 = n1, n2 = n2, p = p,
                    maximum = TRUE, tol = 1e-6 )
  lik2 <- mod2$objective
  k1 <- mod1$maximum   ;   k2 <- mod2$maximum

  while ( abs( lik2 - lik1 ) > tol ) {
    lik1 <- lik2
    mu <- k1 * sy1 + k2 * sy2
    mu <- mu / sqrt( sum(mu^2) )
    R1 <- sum( as.vector(y1 %*% mu) )
    R2 <- sum( as.vector(y2 %*% mu) )

    mod1 <- optimize( vmfa, c(0, 7000), k2 = k2, R1 = R1, R2 = R2, n1 = n1, n2 = n2, p = p,
                      maximum = TRUE, tol = 1e-6 )
    mod2 <- optimize( vmfb, c(0, 7000), k1 = mod1$maximum, R1 = R1, R2 = R2, n1 = n1, n2 = n2, p = p,
                      maximum = TRUE, tol = 1e-6 )
    lik2 <- mod2$objective
    k1 <- mod1$maximum   ;   k2 <- mod2$maximum
  }

  list(mu = mu, kappa1 = k1, kappa2 = k2, loglik = lik2 - 0.5 * n * p * log(2 * pi) )
}
