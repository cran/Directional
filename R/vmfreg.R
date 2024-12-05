vmfreg <- function(y, x, con = TRUE, xnew = NULL, tol = 1e-06) {

  x <- model.matrix( ~., data.frame(x) )
  if ( !con ) x <- x[, -1]
  dm <- dim(y)
  n <- dm[1]  ;  d <- dm[2]

  regvmf <- function(be, y, x, d) {
    be <- matrix( be, ncol = d)
    mu <- x %*% be
    ki <- sqrt( Rfast::rowsums(mu^2) )
    - (0.5 * d - 1) * sum( log(ki) ) - sum(mu * y) +
    sum( log(besselI(ki, d/2 - 1, expon.scaled = TRUE) ) + ki )
  }

  p <- dim(x)[2]

  tic <- proc.time()
  ini <- solve( crossprod(x), crossprod(x, y) )  ## initial values for the beta
  suppressWarnings({
  val1 <- nlm(regvmf, ini, y = y, x = x, d = d,iterlim = 1000)
  val2 <- nlm(regvmf, val1$estimate, y = y, x = x, d = d, iterlim = 1000)
  while (val1$minimum - val2$minimum > tol) {
    val1 <- val2
    val2 <- nlm(regvmf, val1$estimate, y = y, x = x, d = d, iterlim = 1000)
  }
  da <- optim(val2$estimate, regvmf, y = y, x = x, d = d, control = list(maxit = 10000), hessian = TRUE)
  })
  runtime <- proc.time() - tic

  be <- matrix(da$par, ncol = d)
  seb <- sqrt( diag( solve(da$hessian) ) )
  seb <- matrix(seb, ncol = d)

  if ( is.null(xnew) ) {
    mu <- x %*% be
    ki <- sqrt( Rfast::rowsums(mu^2) )
    est <- mu / ki
    fit <- sum( y * est )
  } else {
    xnew <- model.matrix( ~., data.frame(xnew) )
    if ( !con ) xnew <- xnew[, -1]
    mu <- xnew %*% be
    est <- mu / sqrt( Rfast::rowsums(mu^2) )
    fit <- NULL
    ki <- NULL
  }

  if ( is.null( colnames(y) ) ) {
    colnames(est) <- colnames(be) <- colnames(seb) <- paste("Y", 1:d, sep = "")
  } else  colnames(est) <- colnames(be) <- colnames(seb) <- colnames(y)
  rownames(be) <- rownames(seb) <- colnames(x)

  list(runtime = runtime, loglik = -da$value - n * 0.5 * d * log(2 * pi), fit = fit, beta = be, seb = seb, ki = ki, est = est)
}

