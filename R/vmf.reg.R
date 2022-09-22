## Fisher2 regression
vmf.reg <- function(y, x, con = TRUE, xnew = NULL, tol = 1e-06) {

  x <- model.matrix( ~., data.frame(x) )
  if ( !con ) x <- x[, -1]
  n <- dim(y)[1]

  regvmf <- function(be, y, x) {
    be <- matrix( be, ncol = 3)
    mu <- x %*% be
    ki <- sqrt( Rfast::rowsums(mu^2) )
    - sum( log(ki) ) + sum( log(sinh(ki) ) ) - sum(mu * y)
  }
  p <- dim(x)[2]
  ini <- solve(crossprod(x), crossprod(x, y))  ## initial values for the beta
  suppressWarnings({
  val1 <- nlm(regvmf, ini, y = y, x = x, iterlim = 1000)
  val2 <- nlm(regvmf, val1$estimate, y = y, x = x, iterlim = 1000)
  while (val1$minimum - val2$minimum > tol) {
    val1 <- val2
    val2 <- nlm(regvmf, val1$estimate, y = y, x = x, iterlim = 1000)
  }
  da <- optim(val2$estimate, regvmf, y = y, x = x, control = list(maxit = 10000), hessian = TRUE)
  })
  be <- matrix(da$par, ncol = 3)
  seb <- sqrt( diag( solve(da$hessian) ) )
  seb <- matrix(seb, ncol = 3)

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
    colnames(est) <- colnames(be) <- colnames(seb) <- c("X", "Y", "Z")
  } else  colnames(est) <- colnames(be) <- colnames(seb) <- colnames(y)
  rownames(be) <- rownames(seb) <- colnames(x)

  list(loglik = -da$value - n * log(4 * pi), fit = fit, beta = be, seb = seb, ki = ki, est = est)
}
