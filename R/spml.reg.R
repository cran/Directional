spml.reg <- function(y, x, rads = TRUE, xnew = NULL, seb = TRUE, tol = 1e-07) {
  ## y is the angular dependent variable
  ## x contains the independent variable(s)
  ## xnew is some new data or the current ones
  ## pred is either TRUE (xnew is new data) or
  ## FALSE (xnew is the same as x)
  ## if the data are in degrees we transform them into radians
  tic <- proc.time()
  mod <- Rfast::spml.reg(y, x, tol, seb = seb)
  runtime <- proc.time() - tic

  if ( !is.null(xnew) ) {  ## predict new values?
    xnew <-  model.matrix(~., data.frame(xnew) )
    est <- xnew %*% mod$be
  } else {
    x <-  model.matrix(~., data.frame(x) )
    est <- x %*% mod$be
  }
  est <- ( atan(est[, 2]/est[, 1]) + pi * I(est[, 1] < 0) ) %% (2 * pi)
  if ( !rads )  est = est * 180 / pi
  list(runtime = runtime, iters = mod$iter, beta = mod$be, seb = mod$seb, loglik = mod$loglik, est = est)
}
