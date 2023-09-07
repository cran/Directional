sipc.reg <- function(y, x, con = TRUE, xnew = NULL, tol = 1e-06) {
  ## y is the spherical data, unit vectors
  ## x is the independent variables
  ## If con == TRUE a constant term will be estimated estimated
  ini <- Directional::iag.reg(y, x, con = con)$be  ## initial values for the betas

  x <- model.matrix( ~., data.frame(x) )
  if ( !con ) x <- x[, -1]
  n <- dim(y)[1]

  regipc <- function(be, y, x) {
    be <- matrix(be, ncol = 3)
    mu <- x %*% be
    a <- rowSums( y * mu )
    rl <- rowSums(mu^2)
    d <- rl + 1 - a^2
    sqd <- sqrt(d)
    up <- log( (rl + 1) * sqd * ( atan2(sqd, -a) - atan2(sqd, a) + pi ) + 2 * a * d )
    - sum(up) + 2 * sum( log(d) )
  }

  suppressWarnings({
    val1 <- optim( ini, regipc, y = y, x = x, control = list(maxit = 5000) )
    val2 <- optim( val1$par, regipc, y = y, x = x, control = list(maxit = 5000) )
    while ( val1$value - val2$value > tol) {
      val1 <- val2
      val2 <- optim( val1$par, regipc, y = y, x = x, control = list(maxit = 5000) )
    }
    da <- optim(val2$par, regipc, y = y, x = x, control = list(maxit = 10000), hessian = TRUE)
  })

  be <- matrix(da$par, ncol = 3)
  seb <- sqrt( diag( solve(da$hessian) ) )
  seb <- matrix(seb, ncol = 3)

  if ( is.null(xnew) ) {
    mu <- x %*% be
    ki <- sqrt( Rfast::rowsums(mu^2) )
    est <- mu / ki
    fit <- sum( y * est )
    names(fit) <- c( "Fit value")
  } else {
    xnew <- model.matrix( ~., data.frame(xnew) )
    if ( !con ) xnew <- xnew[, -1]
    mu <- xnew %*% be
    est <- mu / sqrt( Rfast::rowsums(mu^2) )
    fit <- NULL
  }

  if ( is.null( colnames(y) ) ) {
    colnames(est) <- colnames(be) <- colnames(seb) <- c("X", "Y", "Z")
  } else  colnames(est) <- colnames(be) <- colnames(seb) <- colnames(y)
  rownames(be) <- rownames(seb) <- colnames(x)

  list(loglik = -val2$value - n * log(4 * pi^2), fit = fit, beta = be, seb = seb, ki = ki, est = est)
}
