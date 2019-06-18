################################
#### MLE estimation of the Kent distribution on the sphere
#### Tsagris Michail 05/2014
#### mtsagris@yahoo.gr
#### References: John Kent(1982)
#### JRSSB, 44(1): 71-80.
#### The Fisher-Bingham distribution on the sphere
################################
kent.mle <- function(x) {
  ## x is the data in Euclidean coordinates
  tic <- proc.time()
  n <- dim(x)[1]  ## sample size
  if ( bigmemory::is.big.matrix(x) ) {
    xbar <- colMeans(x[])
    S <- crossprod(x[]) / n
  } else {
    xbar <- Rfast::colmeans(x)  ## mean vector
    S <- crossprod(x) / n
  }
  xbar <- xbar / sqrt( sum(xbar^2) ) ## mean direction
  u <- c( acos(xbar[1]), ( atan(xbar[3] / xbar[2]) + pi * I(xbar[2]<0) )
  %% (2 * pi) )
  ## u is the mean vector to latitude and longitude
  theta <- u[1]
  phi <- u[2]
  costheta <- cos(theta)
  sintheta <- sin(theta)
  cosphi <- cos(phi)
  sinphi <- sin(phi)
  H <- matrix( c(costheta, sintheta * cosphi, sintheta * sinphi, -sintheta, costheta * cosphi,
       costheta * sinphi, 0, -sinphi, cosphi), ncol = 3 )
  B <- crossprod(H, S) %*% H
  psi <- 0.5 * atan(2 * B[2, 3]/(B[2, 2] - B[3, 3]))
  K <- matrix( c(1, 0, 0, 0, cos(psi), sin(psi), 0, -sin(psi), cos(psi) ), ncol = 3)
  G <- H %*% K  ## The G matrix Kent describes, the A in our notation
  lam <- eigen(B[-1, -1])$values
  ## the next function will be used to estimate the kappa and beta
  if ( bigmemory::is.big.matrix(x) ) {
    options( bigalgebra.mixed_airthmetic_returns_R_matrix = FALSE )
    xg <- x %*% G
    xg1 <- bigmemory::sub.big.matrix(xg, firstCol = 1,lastCol = 1)
    xg1 <- sum(xg1[])
    xg23 <- bigmemory::sub.big.matrix(xg, firstCol = 2,lastCol = 3)
    a <- colSums(xg23[]^2)
    options(bigalgebra.mixed_airthmetic_returns_R_matrix = TRUE)
  } else {
    xg <- x %*% G
    xg1 <- sum(xg[, 1])
    a <- Rfast::colsums(xg[, 2:3]^2)
  }
  xg2 <- a[1]
  xg3 <- a[2]

  mle <- function(para) {
    ## maximization w.r.t. to k and b
    k <- para[1]
    b <- para[2]
    gam <- c(0, k, 0)
    lam <- c(0, -b, b)
    ckb <- fb.saddle(gam, lam)[3]
    g <-  n * ckb - k * xg1 - b * ( xg2 - xg3 )
    g
  }

  ini <- Rfast::vmf.mle(x)$kappa
  ini <- c(ini, ini/2.1)  ## initial values for kappa and beta
  qa <- optim(ini, mle)
  para <- qa$par
  k <- para[1]
  b <- para[2]  ## the estimated parameters
  gam <- c(0, k, 0)
  lam <- c(0, -b, b)
  ckb <- as.numeric( fb.saddle(gam, lam)[3] )
  ## the line below calculates the log-likelihood
  l <-  -n * ckb + k * xg1 + b * ( xg2 - xg3 )
  param <- c(k, b, psi)
  runtime <- proc.time() - tic
  names(param) <- c("kappa", "beta", "psi")
  colnames(G) <- c("mean", "major", "minor")
  list(G = G, param = param, logcon = ckb, loglik = l, runtime = runtime)
}
