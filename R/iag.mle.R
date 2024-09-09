# Function to find MLEs
# @param X, dataset with n observations, n*d matrix
# @return list that contains a vector of MLEs, maxmimum log-likelihood value, etc from function optim() with 'BFGS'. It outputs Hessian matrix as well
iag.mle <- function(x, tol = 1e-6) {

  d <- dim(x)[2]
  if ( d == 3 ) {
    res <- Rfast::iag.mle(x, tol = tol)
  } else if ( d > 3 ) {

    p <- d - 1
    n <- dim(x)[1]

    logCd <- log( (2 * pi)^(-0.5 * p) )
    Mp <- matrix(1, n, p)

    joint_log_lik <- function(mu, x, p, Mp) {
      a <- as.vector( x %*% mu )
      g2 <- sum(mu^2)
      Mp[, 1] <- a * pnorm(a) + dnorm(a)
      Mp[, 2] <- (1 + a^2) * pnorm(a) + a * dnorm(a)
      if ( p >= 3 )  for ( j in 3:p )  Mp[, j] <- a * Mp[, j - 1] + (j - 1) * Mp[, j - 2]
      f <- 0.5 * ( a^2 - g2) + log(Mp[, p])
      - sum(f)
    }
    mod <- optim( par = rep(1, d), joint_log_lik, x = x, p = p, Mp = Mp, method = "BFGS",
                  control = list(maxit = 10000) )

    mesi <- mod$par
    g2 <- sum(mesi^2)
    mesi <- rbind(mesi, mesi / sqrt(g2) )
    rownames(mesi) <- c("Mean vector", "Mean direction")
    param <- c(g2, -mod$value + n * logCd )
    names(param) <- c("Norm of mean", "Log likelihood")
    res <- list( mesi = mesi, param = param )
  }
  res
}
