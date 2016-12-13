################################
#### Projected multivariate normal for circular or angular regression
#### Tsagris Michail 1/2014
#### mtsagris@yahoo.gr
#### Presnell, Morrison and Littell (1998), JASA
################################

spml.reg <- function(y, x, rads = TRUE, xnew = NULL, seb = TRUE) {
  ## y is the angular dependent variable
  ## x contains the independent variable(s)
  ## xnew is some new data or the current ones
  ## pred is either TRUE (xnew is new data) or
  ## FALSE (xnew is the same as x)
  ## if the data are in degrees we transform them into radians
  x <- model.matrix(~., data.frame(x) )
  if ( !rads )   y <- y/180 * pi
  u <- cbind( cos(y), sin(y) )  ## bring the data onto the circle
  n <- nrow(u)
  csx <- crossprod(x)
  XX <- solve( csx, t(x) )
  tXX <- t(XX)
  p <- dim(x)[2]

  funa <- function(be) {
    mu <- x %*% be
    tau <- rowSums(u * mu)
    ell <-  - 0.5 * sum( mu * mu ) + sum( log( 1 + tau * pnorm(tau) / dnorm(tau) ) )
    ell
  }

  tic <- proc.time()
  para <- as.vector( coef( lm.fit(x, u) ) )  ## starting values
  ### E-M algorithm is implemented below
  B <- matrix(para, ncol = 2)
  mu <- x %*% B
  tau <- Rfast::rowsums(u * mu)
  ptau <- pnorm(tau)
  lik1 <-  - 0.5 * sum( mu * mu ) + sum( log( 1 + tau * ptau / dnorm(tau) ) )
  psit <- tau + ptau / ( dnorm(tau) + tau * ptau )
  B <- crossprod( tXX * psit, u)
  mu <- x %*% B
  tau <- Rfast::rowsums(u * mu)
  ptau <- pnorm(tau)
  lik2 <-  - 0.5 * sum( mu * mu ) + sum( log( 1 + tau * ptau / dnorm(tau) ) )

  i <- 2
  while ( lik2 - lik1 > 1e-07 ) {
    i <- i + 1
    lik1 <- lik2
    psit <- tau + ptau / ( dnorm(tau) + tau * ptau )
    B <- crossprod( tXX * psit, u)
    mu <- x %*% B
    tau <- Rfast::rowsums(u * mu)
    ptau <- pnorm(tau)
    lik2 <-  - 0.5 * sum( mu * mu ) + sum( log( 1 + tau * ptau / dnorm(tau) ) )
  }

  loglik <- lik2 - n * log(2 * pi)

  if ( seb ) {

    dtau <- dnorm(tau)
    pdtau <- tau * ptau
    frac <-  ptau/( dtau + pdtau )
    psit <- tau + frac
    psit2 <- 2 - pdtau / (dtau + pdtau)  - frac^2
    C <- u[, 1]    ;    S <- u[, 2]

    A1 <-  - csx
    A2 <-  t( x * psit2 )
    s11 <-  A1 + A2 %*% tcrossprod(C) %*% x
    s12 <- t( tcrossprod(C, A2) ) %*% crossprod(S, x)
    s21 <- t(s12)
    s22 <-  A1 + A2 %*% tcrossprod(S) %*% x
    se1 <- cbind(s11, s12)
    se2 <- cbind(s21, s22)
    se <-  - rbind(se1, se2)  ## negative Hessian of the log-likelihood
    se <- solve(se)
    se <- sqrt( diag(se) )  ## standard errors of the coefficients
    seb <- matrix(se, ncol = 2)
    colnames(seb) <- c("Cosinus of y", "Sinus of y")
    rownames(seb) <- colnames(x)

  } else seb <- NULL

  colnames(B) <- c("Cosinus of y", "Sinus of y")
  rownames(B) <- colnames(x)

  runtime <- proc.time() - tic

  if ( !is.null(xnew) ) {  ## predict new values?
    xnew <-  model.matrix(~., data.frame(xnew) )
    est <- xnew %*% B
    est <- ( atan(est[, 2]/est[, 1]) + pi * I(est[, 1] < 0) ) %% (2 * pi)

  } else  est <- ( atan(mu[, 2]/mu[, 1]) + pi * I(mu[, 1] < 0) ) %% (2 * pi)

  if ( !rads )  est = est * 180 / pi

  list(runtime = runtime, iters = i, beta = B, seb = seb, loglik = loglik, est = est)

}

