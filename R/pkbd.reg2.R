pkbd.reg2 <- function(y, x, con = TRUE, xnew = NULL, tol = 1e-6) {

  x <- model.matrix( ~., data.frame(x) )
  if ( !con )  x <- x[, -1]
  dm <- dim(y)
  n <- dm[1]  ;  d <- dm[2]

  regpkbd <- function(be, y, x, n, d) {
    be <- matrix(be, ncol = d)
    mu <- x %*% be
    a <- Rfast::rowsums(y * mu)
    g2 <- Rfast::rowsums(mu^2)
    com <- sqrt(g2 + 1)
    com2 <- com - a
    rho <- ( com - 1 ) / sqrt(g2)
    0.5 * d * sum( log( com2 ) ) + 0.5 * (d - 2) * sum( log(1 - rho^2) )
  }

  tic <- proc.time()
  ini <- solve(crossprod(x), crossprod(x, y))  ## initial values for the beta
  qa <- optim( as.vector(ini), regpkbd, y = y, x = x, n = n, d = d,
               control = list(maxit = 50000), method = "BFGS" )
  qa <- optim( qa$par, regpkbd, y = y, x = x, n = n, d = d, control = list(maxit = 50000), method = "BFGS" )
  lik1 <- qa$value
  qa <- optim( qa$par, regpkbd, y = y, x = x, n = n, d = d, control = list(maxit = 50000) )
  lik2 <- qa$value

  while ( lik1 - lik2 > tol ) {
    qa <- optim( qa$par, regpkbd, y = y, x = x, n = n, d = d, control = list(maxit = 50000) )
    lik1 <- qa$value
    qa <- optim( qa$par, regpkbd, y = y, x = x, n = n, d = d, control = list(maxit = 50000), method = "BFGS" )
    lik2 <- qa$value
  }

  qa <- optim( qa$par, regpkbd, y = y, x = x, n = n, d = d,
               control = list(maxit = 50000), hessian = TRUE )

  runtime <- proc.time() - tic

  be <- matrix(qa$par, ncol = d)
  seb <- sqrt( diag( solve(qa$hessian) ) )
  seb <- matrix(seb, ncol = d)

  if ( is.null(xnew) ) {
    mu <- x %*% be
    g2 <- sqrt( Rfast::rowsums(mu^2) )
    est <- mu / g2
    fit <- sum( y * est )
  } else {
    xnew <- model.matrix( ~., data.frame(xnew) )
    if ( !con ) xnew <- xnew[, -1]
    mu <- xnew %*% be
    est <- mu / sqrt( Rfast::rowsums(mu^2) )
    fit <- NULL
    g2 <- NULL
  }

  if ( is.null( colnames(y) ) ) {
    colnames(est) <- colnames(be) <- colnames(seb) <- paste("Y", 1:d, sep = "")
  } else  colnames(est) <- colnames(be) <- colnames(seb) <- colnames(y)
  rownames(be) <- rownames(seb) <- colnames(x)

  list( runtime = runtime, loglik = -qa$value + n * lgamma(0.5 * d) - n * 0.5 * d * log(pi) - n * log(2),
        fit = fit, be = be, seb = seb, g2 = g2, est = est)
}
