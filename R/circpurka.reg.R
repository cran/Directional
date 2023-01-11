circpurka.reg <- function(y, x, rads = TRUE, xnew = NULL) {
  if ( !is.matrix(y) ) {
    if ( !rads )   y <- y * pi/180
    z <- cbind( cos(y), sin(y) )
  } else z <- y
  x <- model.matrix( ~., data.frame(x) )

  reg <- function(be, z, x) {
    be <- matrix(be, ncol = 2)
    est <- x %*% be
    a <- sqrt( Rfast::rowsums(est^2) )
    est <- est / a
    - sum( log(a) - log(1 - exp(-a * pi) )  - a * acos( z * est ) )
  }

  ini <- as.vector( solve( crossprod(x), crossprod(x, z) ) )
  mod <- optim(ini, reg, z = z, x = x )
  lik1 <- mod$value
  mod <- optim(mod$par, reg, z = z, x = x, hessian = TRUE )
  lik2 <- mod$value
  while ( lik1 - lik2 > 1e-6 ) {
    lik1 <- lik2
    mod <- optim(mod$par, reg, z = z, x = x, hessian = TRUE )
    lik2 <- mod$value
  }
  be <- matrix(mod$par, ncol = 2)
  covb <- solve( mod$hessian )

  est <- NULL
  if ( !is.null(xnew) ) {
    xnew <- model.matrix(~., data.frame(xnew))
    est <- xnew %*% be
    est <- ( atan(est[, 2]/est[, 1]) + pi * I(est[, 1] < 0) ) %% (2 * pi)
    if ( !rads )  est <- est * 180 / pi
  }
  colnames(be) <- c("Cosinus of y", "Sinus of y")
  rownames(be) <- colnames(x)

  list( value = - mod$value - dim(x)[1] * log(2), be = be, covb = covb, est = est )
}
