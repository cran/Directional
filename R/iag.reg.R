### Regression with the projected normal, no gammas and identity covariance matrix
### Independent angular gaussian normal (rotational symmetry) regression

iag.reg <- function(y, x, con = TRUE, xnew = NULL, tol = 1e-06) {
  ## y is the spherical data, unit vectors
  ## x is the independent variables
  ## If con == TRUE a constant term will be estimated estimated
  x <- model.matrix( ~., data.frame(x) )
  if ( !con ) x <- x[, -1]
  n <- dim(y)[1]

   regiag <- function(be, y, x) {
     be <- matrix(be, ncol = 3)
     mu <- x %*% be
     a <- rowSums( y * mu )
     M2 <- a + (1 + a^2) * pnorm(a)/dnorm(a)
     0.5 * sum(mu * mu) - sum( log(M2) )
   }

  ini <- solve(crossprod(x), crossprod(x, y))  ## initial values for the beta
  oop <- options(warn = -1) 
  val1 <-  nlm(regiag, ini, y = y, x = x, iterlim = 1000)
  val2 <-  nlm(regiag, val1$estimate, y = y, x = x, iterlim = 1000)
  while (val1$minimum - val2$minimum > tol) {
    val1 <- val2
    val2 <-  nlm(regiag, val1$estimate, y = y, x = x, iterlim = 1000)
  }
  da <- optim(val2$estimate, regiag, y = y, x = x, control = list(maxit = 10000), hessian = TRUE)
  on.exit( options(oop) )
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
  }

  if ( is.null( colnames(y) ) ) {
    colnames(est) <- colnames(be) <- colnames(seb) <- c("X", "Y", "Z")
  } else  colnames(est) <- colnames(be) <- colnames(seb) <- colnames(y)
  rownames(be) <- rownames(seb) <- colnames(x)
  names(fit) <- c( "Fit value")

  list(loglik = -da$value - 1.5 * n * log(2 * pi), fit = fit, beta = be, seb = seb, ki = ki, est = est)
}

