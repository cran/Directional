circpurka.reg <- function(y, x, rads = TRUE, xnew = NULL, tol = 1e-6, maxiters = 300) {

  tic <- proc.time()
  if ( !is.matrix(y) ) {
    if ( !rads )   y <- y * pi/180
    z <- cbind( cos(y), sin(y) )
  } else z <- y
  x <- model.matrix( ~., data.frame(x) )

  suppressWarnings({
    ini <- .reg.nr(z, x, tol = tol, maxiters = maxiters)
    mod1 <- optim(ini, .reg, z = z, x = x )
    lik1 <- mod1$value
    mod2 <- try( optim(mod1$par, .reg, z = z, x = x, hessian = TRUE ), silent = TRUE )
    if ( identical( class(mod2), "try-error" ) ) {
      mod2 <- mod1
    } else  lik2 <- mod2$value
    while ( mod1$value - mod2$value > 1e-6 ) {
      mod1 <- mod2
      mod2 <- try( optim(mod1$par, .reg, z = z, x = x, hessian = TRUE ), silent = TRUE )
      if ( identical( class(mod2), "try-error" ) )  mod2 <- mod1
    }
  })
  be <- matrix(mod2$par, ncol = 2)
  colnames(be) <- c("Cosinus of y", "Sinus of y")
  rownames(be) <- colnames(x)
  seb <- NULL
  if ( !is.null(mod2$hessian) ) {
    seb <- solve(mod2$hessian)
    seb <- matrix( sqrt( diag(seb) ), ncol = 2 )
    colnames(seb) <- c("Cosinus of y", "Sinus of y")
    rownames(seb) <- colnames(x)
  }

  runtime <- proc.time() - tic

  est <- NULL
  if ( !is.null(xnew) ) {
    xnew <- model.matrix(~., data.frame(xnew))
    est <- xnew %*% be
    est <- ( atan(est[, 2]/est[, 1]) + pi * I(est[, 1] < 0) ) %% (2 * pi)
    if ( !rads )  est <- est * 180 / pi
  }

  list( runtime = runtime, be = be, seb = seb, loglik = - mod2$value - dim(x)[1] * log(2), est = est )
}




.reg <- function(be, z, x) {
  be <- matrix(be, ncol = 2)
  est <- x %*% be
  a <- sqrt( Rfast::rowsums(est^2) )
  est <- est / a
  d <- Rfast::rowsums(z * est)
  d <- pmin( pmax(d, -1 + 1e-20), 1 - 1e-20 )
  - sum( log(a) - log(1 - exp(-pi * a) ) - a * acos( d ) )
}

.reg.score <- function(be, z, x) {
  be  <- matrix(be, ncol = 2)
  eta <- x %*% be
  a <- sqrt( Rfast::rowsums(eta^2) )
  m <- eta / a
  c <- Rfast::rowsums(z * m)
  c <- pmin( pmax(c, -1 + 1e-20), 1 - 1e-20 )
  s <- sqrt(1 - c^2)
  d <- acos(c)
  phi <- pi / expm1(pi * a)
  w <- 1 / a - phi - d
  tang <- (z - c * m) / s
  G <- w * m + tang                     # n x 2, grad of loglik_i wrt eta_i
  cbind(-G[, 1] * x, -G[, 2] * x)          # n x 2p, grad of -loglik_i wrt be
}

.reg.nr <- function(z, x, tol = 1e-6, maxiters = 300) {

  be <- as.vector( solve( crossprod(x), crossprod(x, z) ) )

  f <- function(b)  .reg(b, z, x)
  l <- f(be)

  for ( i in 1:maxiters ) {
    com <- .reg.score(be, z, x)
    g <-  Rfast::colsums(com)
    I <- crossprod(com)
    step <- tryCatch( solve(I, g), error = function(e) g )  # fallback if singular
    a.step <- 1
    be.new <- be - a.step * step
    l.new  <- f(be.new)
    while ( !is.finite(l.new) || l.new > l ) {
      a.step <- a.step / 2
      if (a.step < 1e-14) { be.new <- be; l.new <- l; break }
      be.new <- be - a.step * step
      l.new  <- f(be.new)
    }
    conv <- abs(l.new - l) < tol
    be <- be.new ;  l <- l.new
    if (conv)  break
  }
  be
}

