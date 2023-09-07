cipc.reg <- function(y, x, rads = TRUE, xnew = NULL) {

  lik <- function(be, y, x) {
    be <- matrix(be, ncol = 2)
    mu <- x %*% be
    g2 <- rowsums(mu^2)
    a <- rowsums(y * mu)
    sum( log(sqrt(g2 + 1) - a) )
  }

  if ( !is.matrix(y) ) {
    if ( !rads )  y <- y * pi / 180
    y <- cbind( cos(y), sin(y) )
  }

  x <- model.matrix(y~., as.data.frame(x), )
  n <- dim(y)[1]  ;  p <- dim(x)[2]
  runtime <- proc.time()
  mod <- optim( rnorm(2 * p), lik, y = y, x = x, control = list(maxit = 5000) )
  lik1 <-  - mod$value
  mod <- optim( mod$par, lik, y = y, x = x, control = list(maxit = 5000) )
  lik2 <-  - mod$value
  while (lik2 - lik1 > 1e-6) {
    lik1 <- lik2
    mod <- optim( mod$par, lik, y = y, x = x, control = list(maxit = 5000) )
    lik2 <-  -mod$value
  }
  mod <- optim( mod$par, lik, y = y, x = x, control = list(maxit = 5000),
                hessian = TRUE )

  runtime <- proc.time() - runtime
  be <- mod$par
  be <- matrix(be, ncol = 2)
  seb <- solve( mod$hessian )
  seb <- matrix( sqrt( diag(seb) ), ncol = 2)

  colnames(be) <- c("Cosinus of y", "Sinus of y")
  rownames(be) <- colnames(x)
  colnames(seb) <- c("Cosinus of y", "Sinus of y")
  rownames(seb) <- colnames(x)

  est <- NULL
  if ( !is.null(xnew) ) {
    xnew <- model.matrix(~., data.frame(xnew))
    est <- xnew %*% be
    est <- ( atan(est[, 2]/est[, 1]) + pi * I(est[, 1] < 0) ) %% (2 * pi)
    if ( !rads )  est <- est * 180/pi
  }

  list( runtime = runtime, beta = be, seb = seb, loglik = -mod$value - n * log(2 * pi),
        est = est )
}
