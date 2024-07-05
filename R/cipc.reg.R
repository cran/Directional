cipc.reg <- function(y, x, rads = TRUE, xnew = NULL, tol = 1e-6) {

  if ( !is.matrix(y) ) {
    if ( !rads )  y <- y * pi / 180
    y <- cbind( cos(y), sin(y) )
  }

  x <- model.matrix(y~., as.data.frame(x) )
  n <- dim(y)[1]  ;  p <- dim(x)[2]
  H <- matrix(0, 2 * p, 2 * p)

  tic <- proc.time()
  be <- Rfast::spml.reg(y, x[, -1], tol)$be
  mu <- x %*% be
  g2 <- Rfast::rowsums(mu^2)
  a <- Rfast::rowsums(y * mu)
  com <- sqrt(g2 + 1)
  com2 <- com - a
  lik <-  - sum( log( com2 ) )

  muc_y <- mu / com - y
  der1 <- Rfast::eachcol.apply(x, muc_y[, 1] / com2 )
  der2 <- Rfast::eachcol.apply(x, muc_y[, 2] / com2 )

  ### Jacobian of b1
  a1 <- ( com - mu[, 1]^2 / com ) / ( com^2 * com2 )
  up1 <- crossprod(x, x * a1)
  up2 <- crossprod(x * muc_y[, 1]/com2)
  H[1:p, 1:p] <- up2 - up1
  ### Jacobian of b2
  a1 <- ( com - mu[, 2]^2 / com ) / ( com^2 * com2 )
  up1 <- crossprod(x, x * a1)
  up2 <- crossprod(x * muc_y[, 2]/com2)
  H[(p + 1):(2*p), (p + 1):(2*p)] <- up2 - up1
  ### Jacobian of b12
  a1 <- mu[, 1] * mu[, 2] / ( com^3 * com2)
  up1 <- crossprod(x, x * a1)
  up2 <- crossprod(x * muc_y[, 1]/com2, x * muc_y[, 2]/com2)
  H[1:p, (p + 1):(2*p)] <- H[(p + 1):(2*p), 1:p] <- up2 + up1

  be <- be + solve(H, c(der1, der2))
  mu <- x %*% be
  g2 <- Rfast::rowsums(mu^2)
  a <- Rfast::rowsums(y * mu)
  com <- sqrt(g2 + 1)
  com2 <- com - a
  lik[2] <-  - sum( log( com2 ) )

  i <- 2
  while (lik[i] - lik[i-1] > tol) {
    i <- i + 1
    muc_y <- mu / com - y
    der1 <- Rfast::eachcol.apply(x, muc_y[, 1] / com2 )
    der2 <- Rfast::eachcol.apply(x, muc_y[, 2] / com2 )

    ### Jacobian of b1
    a1 <- ( com - mu[, 1]^2 / com ) / ( com^2 * com2 )
    up1 <- crossprod(x, x * a1)
    up2 <- crossprod(x * muc_y[, 1]/com2)
    H[1:p, 1:p] <- up2 - up1
    ### Jacobian of b2
    a1 <- ( com - mu[, 2]^2 / com ) / ( com^2 * com2 )
    up1 <- crossprod(x, x * a1)
    up2 <- crossprod(x * muc_y[, 2]/com2)
    H[(p + 1):(2*p), (p + 1):(2*p)] <- up2 - up1
    ### Jacobian of b12
    a1 <- mu[, 1] * mu[, 2] / ( com^3 * com2)
    up1 <- crossprod(x, x * a1)
    up2 <- crossprod(x * muc_y[, 1]/com2, x * muc_y[, 2]/com2)
    H[1:p, (p + 1):(2*p)] <- H[(p + 1):(2*p), 1:p] <- up2 + up1

    be <- be + solve(H, c(der1, der2))
    mu <- x %*% be
    g2 <- Rfast::rowsums(mu^2)
    a <- Rfast::rowsums(y * mu)
    com <- sqrt(g2 + 1)
    com2 <- com - a
    lik[i] <-  - sum( log( com2 ) )
  }

  runtime <- proc.time() - tic
  seb <- solve( - H )
  seb <- matrix( sqrt( diag(seb) ), ncol = 2)
  colnames(be) <- colnames(seb) <- c("Cosinus of y", "Sinus of y")
  rownames(be) <- rownames(seb) <- colnames(x)

  est <- NULL
  if ( !is.null(xnew) ) {
    xnew <- model.matrix(~., data.frame(xnew))
    est <- xnew %*% be
    est <- ( atan(est[, 2]/est[, 1]) + pi * I(est[, 1] < 0) ) %% (2 * pi)
    if ( !rads )  est <- est * 180/pi
  }

  list( runtime = runtime, iters = i, beta = be, seb = seb,
        loglik = lik[i] - n * log(2 * pi), est = est )
}









# cipc.reg <- function(y, x, rads = TRUE, xnew = NULL) {
#
#   lik <- function(be, y, x) {
#     be <- matrix(be, ncol = 2)
#     mu <- x %*% be
#     g2 <- Rfast::rowsums(mu^2)
#     a <- Rfast::rowsums(y * mu)
#     sum( log(sqrt(g2 + 1) - a) )
#   }
#
#   if ( !is.matrix(y) ) {
#     if ( !rads )  y <- y * pi / 180
#     y <- cbind( cos(y), sin(y) )
#   }
#
#   x <- model.matrix(y~., as.data.frame(x), )
#   n <- dim(y)[1]  ;  p <- dim(x)[2]
#   runtime <- proc.time()
#   mod <- optim( rnorm(2 * p), lik, y = y, x = x, control = list(maxit = 5000) )
#   lik1 <-  - mod$value
#   mod <- optim( mod$par, lik, y = y, x = x, control = list(maxit = 5000) )
#   lik2 <-  - mod$value
#   while (lik2 - lik1 > 1e-6) {
#     lik1 <- lik2
#     mod <- optim( mod$par, lik, y = y, x = x, control = list(maxit = 5000) )
#     lik2 <-  -mod$value
#   }
#   mod <- optim( mod$par, lik, y = y, x = x, control = list(maxit = 5000),
#                 hessian = TRUE )
#
#   runtime <- proc.time() - runtime
#   be <- mod$par
#   be <- matrix(be, ncol = 2)
#   seb <- solve( mod$hessian )
#   seb <- matrix( sqrt( diag(seb) ), ncol = 2)
#
#   colnames(be) <- colnames(seb) <- c("Cosinus of y", "Sinus of y")
#   rownames(be) <- rownames(seb) <- colnames(x)
#
#   est <- NULL
#   if ( !is.null(xnew) ) {
#     xnew <- model.matrix(~., data.frame(xnew))
#     est <- xnew %*% be
#     est <- ( atan(est[, 2]/est[, 1]) + pi * I(est[, 1] < 0) ) %% (2 * pi)
#     if ( !rads )  est <- est * 180/pi
#   }
#
#   list( runtime = runtime, beta = be, seb = seb, loglik = -mod$value - n * log(2 * pi),
#         est = est )
# }
